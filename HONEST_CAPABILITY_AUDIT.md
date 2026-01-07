# biometal: Honest Capability Audit
**Generated**: January 7, 2026
**Purpose**: Document what actually works vs what's advertised

## ✅ **WHAT ACTUALLY WORKS** (Tested & Verified)

### File Format Processing
- ✅ **Local FASTQ/FASTA files** - Compressed (.gz) and uncompressed
- ✅ **Format detection** - Auto-detects FASTQ vs FASTA
- ✅ **Large file handling** - Tested on 340MB+ files
- ✅ **Gzip decompression** - Transparent handling

### Core Analysis Commands
- ✅ **count-reads** - Fast, accurate read counting
- ✅ **gc-content** - GC percentage calculation
- ✅ **count-bases** - Base frequency analysis
- ✅ **mean-quality** - Quality score statistics
- ✅ **complexity-score** - Shannon entropy calculation
- ✅ **fastq-to-fasta** - Format conversion

### ARM Performance
- ✅ **NEON acceleration** - 15-25× speedup confirmed on real data
- ✅ **Memory efficiency** - Processes 340MB files with constant memory
- ✅ **Professional output** - Clean formatting, proper statistics

## ❌ **WHAT'S BROKEN** (Advertised but Non-Functional)

### Network/Remote I/O
- ❌ **HTTP/HTTPS URLs** - Returns "No such file or directory"
- ❌ **Network streaming** - Not implemented despite CLI claims
- ❌ **SRA accessions** - Treated as local file paths
- ❌ **stdin support** - Recognizes stdin but processes 0 reads

### Advanced Features
- ❌ **Host removal** - No working alignment against references
- ❌ **Alignment generation** - Statistical scoring integration broke basic tests
- ❌ **Metagenomic workflows** - No reference genome handling

## ⚠️ **PARTIALLY WORKING** (Implemented but Unreliable)

### Statistical Analysis
- ⚠️ **Statistical scoring** - Functions exist but broke integration tests
- ⚠️ **Alignment primitive** - Proof-of-concept exists but not production-ready

## 🎯 **ACTUAL VALUE PROPOSITION**

Based on testing, biometal's **real strengths** are:

1. **Local file processing** with excellent ARM performance
2. **Memory-efficient analysis** of large genomic files
3. **Professional CLI experience** with clean output
4. **Format conversion** and basic statistics

## 📊 **PERFORMANCE VALIDATION**

**Real-world test**: 1.6M read IBD metagenomic dataset (340MB)
- **Base counting**: 16.7× NEON speedup
- **GC content**: 20.3× NEON speedup
- **Quality analysis**: 25.1× NEON speedup
- **Memory usage**: Constant regardless of file size

## ⚠️ **USER RECOMMENDATIONS**

### Use biometal FOR:
- ✅ Processing local FASTQ/FASTA files
- ✅ Fast sequence statistics on ARM hardware
- ✅ Format conversion tasks
- ✅ Memory-constrained environments

### DON'T use biometal FOR:
- ❌ Downloading data from URLs
- ❌ Host removal/contamination filtering
- ❌ Read alignment generation
- ❌ Network-based workflows

## 🔧 **CREDIBILITY RECOMMENDATIONS**

1. **Update CLI help** - Remove claims about URLs, SRA, network streaming
2. **Fix documentation** - Be honest about what works
3. **Focus development** - Double down on proven strengths
4. **Test before claiming** - Don't advertise untested features

## 📈 **CONCLUSION**

biometal **excels** at its core mission: fast, memory-efficient processing of local genomic files on ARM hardware. The ARM NEON acceleration is real and impressive.

The **credibility problem** comes from overpromising features that don't work. Focus on the excellent foundation rather than broken advanced features.

**Grade: B+ for what works, F for truth in advertising**