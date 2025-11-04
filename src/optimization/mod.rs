//! Platform detection and auto-optimization
//!
//! This module implements automatic platform detection and evidence-based
//! threshold selection from OPTIMIZATION_RULES.md.

/// Evidence-based block size for streaming (Rule 2)
///
/// From Entry 027: 10K records balances SIMD efficiency and memory footprint
pub const BLOCK_SIZE: usize = 10_000;

/// Memory-mapped I/O threshold (Rule 4)
///
/// From Entry 032: Use mmap for files ≥50 MB, standard I/O for smaller files
pub const MMAP_THRESHOLD: u64 = 50 * 1024 * 1024; // 50 MB

/// Detect if ARM NEON is available
pub fn has_neon() -> bool {
    cfg!(target_arch = "aarch64")
}

/// Detect platform for optimization selection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Platform {
    /// macOS with Apple Silicon
    MacOS,
    /// Linux on ARM (Graviton, Ampere, etc.)
    LinuxARM,
    /// Linux on x86_64
    LinuxX86,
    /// Windows on ARM
    WindowsARM,
    /// Windows on x86_64
    WindowsX86,
    /// Other/unknown
    Other,
}

impl Platform {
    /// Detect current platform
    pub fn detect() -> Self {
        #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
        return Platform::MacOS;

        #[cfg(all(target_os = "linux", target_arch = "aarch64"))]
        return Platform::LinuxARM;

        #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
        return Platform::LinuxX86;

        #[cfg(all(target_os = "windows", target_arch = "aarch64"))]
        return Platform::WindowsARM;

        #[cfg(all(target_os = "windows", target_arch = "x86_64"))]
        return Platform::WindowsX86;

        #[cfg(not(any(
            all(target_os = "macos", target_arch = "aarch64"),
            all(target_os = "linux", target_arch = "aarch64"),
            all(target_os = "linux", target_arch = "x86_64"),
            all(target_os = "windows", target_arch = "aarch64"),
            all(target_os = "windows", target_arch = "x86_64"),
        )))]
        Platform::Other
    }

    /// Check if mmap optimization is available (Rule 4)
    ///
    /// From Entry 032: Currently validated on macOS only
    pub fn supports_mmap_optimization(self) -> bool {
        matches!(self, Platform::MacOS)
        // Future: Add Platform::LinuxARM after Graviton validation
    }
}
