# Performance Tuning Guide

This guide explains how to optimize biometal's network streaming performance for your specific use case.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Background Prefetching](#background-prefetching)
3. [Cache Configuration](#cache-configuration)
4. [Chunk Size Optimization](#chunk-size-optimization)
5. [Network Characteristics](#network-characteristics)
6. [Memory vs Performance Tradeoffs](#memory-vs-performance-tradeoffs)
7. [Troubleshooting](#troubleshooting)

---

## Quick Start

For most users, the default configuration works well:

```rust
use biometal::io::DataSource;
use biometal::FastqStream;

// Default configuration (recommended for most use cases)
let source = DataSource::Sra("SRR390728".to_string());
let stream = FastqStream::new(source)?;

// Process records...
for record in stream {
    let record = record?;
    // Your analysis here
}
```

**Default settings:**
- Prefetch count: 4 blocks
- Cache size: 50 MB (byte-bounded LRU)
- Chunk size: 64 KB
- Retry attempts: 3 (with exponential backoff)

---

## Background Prefetching

### What is Prefetching?

Background prefetching hides network latency by fetching N blocks ahead in separate threads while you process the current block. When you need the next block, it's already in cache.

### Evidence

**Entry 028 (Lab Notebook):** I/O bottleneck is 264-352× slower than compute. Prefetching is essential for maintaining analysis throughput over networks.

### Configuration

```rust
use biometal::io::{HttpReader, sra_to_url};

let url = sra_to_url("SRR390728")?;
let reader = HttpReader::new(&url)?
    .with_prefetch_count(8);  // Prefetch 8 blocks ahead
```

### Tuning Guidelines

| Network Type | Latency | Bandwidth | Recommended Prefetch |
|-------------|---------|-----------|---------------------|
| Local network / Data center | <5ms | >100 Mbps | 2-4 |
| Home broadband | 10-50ms | 10-100 Mbps | 4-8 (default) |
| Mobile (4G/5G) | 20-100ms | 1-50 Mbps | 8-16 |
| Satellite / 3G | 100-500ms | <1 Mbps | 16-32 |

### Memory Impact

Each prefetch block consumes ~64 KB (default chunk size):
- Prefetch count 2: ~128 KB additional memory
- Prefetch count 4: ~256 KB additional memory
- Prefetch count 8: ~512 KB additional memory
- Prefetch count 16: ~1 MB additional memory

**Note:** This is in addition to the 50 MB cache and ~5 MB streaming buffer.

---

## Cache Configuration

### LRU Cache

biometal uses a byte-bounded LRU (Least Recently Used) cache that automatically evicts old blocks when the memory limit is reached.

### Default: 50 MB

The default 50 MB cache size is chosen to balance:
- **Memory footprint:** Reasonable for most systems
- **Cache hit rate:** Good for sequential and nearby access patterns
- **Eviction overhead:** Minimal for typical workloads

### When to Adjust Cache Size

**Increase cache size (100-200 MB) when:**
- You have abundant RAM (>16 GB)
- Repeated passes over the same data
- Random access patterns
- Working with multiple concurrent streams

**Decrease cache size (10-25 MB) when:**
- Limited RAM (e.g., Raspberry Pi with 4 GB)
- Memory-constrained environments
- Strictly sequential access (prefetching is more important)

### Future API (Planned)

```rust
// Cache size configuration (future API)
let client = HttpClient::new()
    .with_cache_size(200 * 1024 * 1024);  // 200 MB cache
```

---

## Chunk Size Optimization

### Default: 64 KB

The 64 KB chunk size balances:
- **HTTP overhead:** Fewer requests (vs smaller chunks)
- **Latency:** Lower per-request latency (vs larger chunks)
- **Memory:** Reasonable memory per block

### Benchmarking Different Chunk Sizes

Run the network streaming benchmarks to test chunk sizes:

```bash
cargo bench --bench network_streaming --features network
```

The `bench_chunk_size_impact` benchmark tests: 4 KB, 16 KB, 64 KB, 256 KB, 1 MB

### Tuning Guidelines

| Chunk Size | Use Case | Pros | Cons |
|-----------|----------|------|------|
| 4-16 KB | Very high latency networks | Lower per-request latency | More HTTP overhead |
| 64 KB (default) | Balanced for most use cases | Good balance | N/A |
| 256 KB - 1 MB | High bandwidth, low latency | Fewer requests | Higher memory, latency per request |

### Configuration

```rust
let reader = HttpReader::new(&url)?
    .with_chunk_size(128 * 1024);  // 128 KB chunks
```

---

## Network Characteristics

### Measuring Your Network

Use the prefetch tuning example to measure performance:

```bash
cargo run --example prefetch_tuning --features network
```

This tests different prefetch configurations with real SRA data and reports:
- Throughput (reads/sec)
- Memory footprint
- Recommendations based on results

### High Latency Networks

**Symptoms:**
- Slow throughput despite good bandwidth
- Long pauses between processing batches

**Solutions:**
- Increase prefetch count (8-16)
- Increase chunk size (128-256 KB)
- Enable more aggressive caching

### Low Bandwidth Networks

**Symptoms:**
- Consistently slow throughput
- Network saturation

**Solutions:**
- Smaller chunk size (32-64 KB) for lower latency
- Moderate prefetch (4-8) to avoid overwhelming network
- Ensure cache is enabled (50+ MB)

### Unstable Networks

**Symptoms:**
- Intermittent failures
- Variable performance

**Solutions:**
- Default retry logic handles most cases (3 attempts, exponential backoff)
- Increase retry timeout (future API)
- Consider downloading dataset locally if persistently unstable

---

## Memory vs Performance Tradeoffs

### Memory Budget

Calculate your total memory footprint:

```
Total Memory = Streaming Buffer + Cache + Prefetch
             = 5 MB + cache_size + (prefetch_count × chunk_size)
```

**Examples:**

| Configuration | Streaming | Cache | Prefetch | Total |
|--------------|-----------|-------|----------|-------|
| Minimal | 5 MB | 10 MB | 2 × 64 KB | ~15 MB |
| Default | 5 MB | 50 MB | 4 × 64 KB | ~55 MB |
| Aggressive | 5 MB | 200 MB | 16 × 256 KB | ~209 MB |

### Performance Profile

**Best case scenario** (high bandwidth, low latency, sequential access):
- Configuration: Prefetch 4, cache 50 MB, chunk 64 KB
- Achieves near-local performance
- Latency effectively hidden by prefetching

**Worst case scenario** (low bandwidth, high latency, random access):
- Configuration: Prefetch 16, cache 200 MB, chunk 128 KB
- Still maintains constant memory (vs downloading entire dataset)
- Slower than local, but enables analysis of datasets larger than disk/RAM

---

## Troubleshooting

### Slow Performance

**Diagnostic steps:**

1. **Check network connection:**
   ```bash
   # Test connection to NCBI SRA
   curl -I https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR000/SRR000001/SRR000001
   ```

2. **Run benchmarks:**
   ```bash
   cargo bench --bench network_streaming --features network
   ```

3. **Test with real data:**
   ```bash
   cargo run --example prefetch_tuning --features network
   ```

**Common causes:**
- Network latency too high → Increase prefetch count
- Network bandwidth too low → Consider local download
- Cache too small → Increase cache size
- Chunk size suboptimal → Run chunk size benchmark

### High Memory Usage

**Diagnostic steps:**

1. **Check configuration:**
   - What is your cache size? (default 50 MB)
   - What is your prefetch count? (default 4)
   - Are you accumulating records? (should use streaming iterator)

2. **Reduce memory footprint:**
   ```rust
   // Minimal configuration (~15 MB total)
   let reader = HttpReader::new(&url)?
       .with_prefetch_count(2)
       .with_chunk_size(32 * 1024);

   // Note: Reduce cache size when API available
   ```

3. **Verify streaming usage:**
   ```rust
   // ❌ BAD: Accumulates in memory
   let records: Vec<_> = stream.collect();

   // ✅ GOOD: Constant memory
   for record in stream {
       process(record?);
   }
   ```

### Connection Failures

**Built-in retry logic:**
- Automatic retry: 3 attempts
- Exponential backoff: 1s, 2s, 4s
- Handles transient network failures

**Persistent failures:**
- Check SRA accession exists
- Verify internet connection
- Check firewall rules
- Try with different dataset

---

## Benchmarking Your Configuration

### 1. Run Built-in Benchmarks

```bash
# Network streaming benchmarks
cargo bench --bench network_streaming --features network

# Operations benchmarks (ARM NEON)
cargo bench --bench operations
```

### 2. Test with Real Data

```bash
# E. coli dataset (~40 MB)
cargo run --example sra_ecoli --features network

# Test prefetch configurations
cargo run --example prefetch_tuning --features network
```

### 3. Create Custom Benchmarks

```rust
use std::time::Instant;
use biometal::io::DataSource;
use biometal::FastqStream;

fn benchmark_my_analysis() -> biometal::Result<()> {
    let start = Instant::now();

    let source = DataSource::Sra("SRR390728".to_string());
    let stream = FastqStream::new(source)?;

    let mut count = 0;
    for record in stream {
        let record = record?;
        // Your analysis here
        count += 1;
    }

    let elapsed = start.elapsed();
    println!("Processed {} records in {:.2} sec", count, elapsed.as_secs_f64());
    println!("Throughput: {:.1} reads/sec", count as f64 / elapsed.as_secs_f64());

    Ok(())
}
```

---

## Evidence Base

All tuning recommendations are based on experimental validation:

**Entry 028 (Lab Notebook):**
- I/O bottleneck: 264-352× slower than compute
- Network streaming addresses critical bottleneck
- Prefetching + caching essential for performance

**Entry 026:**
- Streaming achieves 99.5% memory reduction
- Constant ~5 MB regardless of dataset size

**See:** [OPTIMIZATION_RULES.md](../OPTIMIZATION_RULES.md) for complete evidence base.

---

## Quick Reference

| Scenario | Prefetch | Cache | Chunk Size |
|----------|----------|-------|-----------|
| Default (recommended) | 4 | 50 MB | 64 KB |
| High latency network | 8-16 | 50-100 MB | 128 KB |
| Low bandwidth network | 4-8 | 50 MB | 32-64 KB |
| Memory constrained | 2 | 10 MB | 32 KB |
| High memory available | 16 | 200 MB | 128 KB |
| Random access | 4 | 100-200 MB | 64 KB |
| Sequential streaming | 8 | 50 MB | 64-128 KB |

---

**Last Updated:** November 4, 2025
**Version:** v0.2.2 (Week 3-4 Polish)
