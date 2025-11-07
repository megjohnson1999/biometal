//! Sliding window minimum using monotonic deque (O(1) amortized)
//!
//! This module implements the "two stacks" or monotonic deque algorithm for finding
//! the minimum value in a sliding window. This is a critical component for efficient
//! minimizer extraction.
//!
//! # Algorithm
//!
//! The algorithm maintains a double-ended queue (deque) with the following properties:
//! 1. Elements are stored with their (value, position) pairs
//! 2. The deque is monotonically non-decreasing in both value and position
//! 3. The front element is always the minimum in the current window
//!
//! When processing a new element:
//! 1. Remove larger elements from the back (they'll never be the minimum)
//! 2. Add the new element to the back
//! 3. Remove elements from the front that fall outside the window
//!
//! # Complexity
//!
//! - **Time**: O(1) amortized per element (each element pushed/popped at most once)
//! - **Space**: O(w) where w is the window size
//!
//! # Evidence
//!
//! - SimdMinimizers: 820 Mbp/s (221× faster than Entry 036 baseline)
//! - Port from simd-minimizers crate (MIT licensed)
//! - Target: 100-200× speedup for minimizers in biometal v1.3.0
//!
//! # Example
//!
//! ```
//! use biometal::operations::sliding_min::{SlidingMin, MinElem};
//!
//! let values = vec![3, 1, 4, 1, 5, 9, 2, 6];
//! let w = 3; // window size
//!
//! let mut sliding_min = SlidingMin::new(w).unwrap();
//! let minima: Vec<MinElem<i32>> = values.into_iter()
//!     .enumerate()
//!     .filter_map(|(pos, val)| sliding_min.push(val, pos))
//!     .collect();
//!
//! // First window [3,1,4]: minimum is 1 at position 1
//! assert_eq!(minima[0].val, 1);
//! assert_eq!(minima[0].pos, 1);
//! ```

use crate::error::{BiometalError, Result};
use std::collections::VecDeque;

/// A value with its absolute position in the input sequence.
///
/// When comparing, ties between values are broken in favor of smaller position
/// (leftmost minimum is preferred).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct MinElem<V> {
    /// The value
    pub val: V,
    /// The absolute position in the input sequence
    pub pos: usize,
}

impl<V: Ord> PartialOrd for MinElem<V> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<V: Ord> Ord for MinElem<V> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Compare by value first, then by position (tie-breaker)
        match self.val.cmp(&other.val) {
            std::cmp::Ordering::Equal => self.pos.cmp(&other.pos),
            other => other,
        }
    }
}

/// Sliding window minimum tracker using monotonic deque.
///
/// This structure efficiently tracks the minimum value in a sliding window
/// of size `w` as new elements are added.
///
/// # Complexity
///
/// - `push`: O(1) amortized (each element pushed/popped at most once)
/// - Space: O(w) worst case
///
/// # Example
///
/// ```
/// use biometal::operations::sliding_min::SlidingMin;
///
/// let mut tracker = SlidingMin::new(3).unwrap();
///
/// // Window: [5]
/// assert_eq!(tracker.push(5, 0), None); // Window not full yet
///
/// // Window: [5, 3]
/// assert_eq!(tracker.push(3, 1), None); // Window not full yet
///
/// // Window: [5, 3, 7]
/// let min1 = tracker.push(7, 2).unwrap();
/// assert_eq!(min1.val, 3); // Minimum of [5, 3, 7]
/// assert_eq!(min1.pos, 1);
///
/// // Window: [3, 7, 2]
/// let min2 = tracker.push(2, 3).unwrap();
/// assert_eq!(min2.val, 2); // Minimum of [3, 7, 2]
/// assert_eq!(min2.pos, 3);
/// ```
pub struct SlidingMin<V> {
    /// Window size
    w: usize,
    /// Monotonic deque: non-decreasing in both val and pos
    /// Front element is always the minimum in the current window
    deque: VecDeque<MinElem<V>>,
    /// Current position (number of elements seen so far)
    current_pos: usize,
}

impl<V: Ord + Copy> SlidingMin<V> {
    /// Create a new sliding window minimum tracker.
    ///
    /// # Arguments
    ///
    /// * `w` - Window size (must be > 0)
    ///
    /// # Errors
    ///
    /// Returns `InvalidInput` if `w` is 0.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::operations::sliding_min::SlidingMin;
    ///
    /// let tracker = SlidingMin::<u64>::new(10).unwrap();
    /// ```
    pub fn new(w: usize) -> Result<Self> {
        if w == 0 {
            return Err(BiometalError::InvalidInput {
                msg: "Window size must be greater than 0".to_string(),
            });
        }

        Ok(SlidingMin {
            w,
            deque: VecDeque::with_capacity(w),
            current_pos: 0,
        })
    }

    /// Push a new value and return the minimum of the current window (if window is full).
    ///
    /// # Arguments
    ///
    /// * `val` - The value to add
    /// * `pos` - The absolute position of this value in the input sequence
    ///
    /// # Returns
    ///
    /// - `Some(MinElem)` if the window is full (returns the minimum)
    /// - `None` if the window is not yet full (first w-1 elements)
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::operations::sliding_min::SlidingMin;
    ///
    /// let mut tracker = SlidingMin::new(3).unwrap();
    ///
    /// assert_eq!(tracker.push(5, 0), None); // Window: [5] (not full)
    /// assert_eq!(tracker.push(3, 1), None); // Window: [5, 3] (not full)
    ///
    /// let min = tracker.push(7, 2).unwrap(); // Window: [5, 3, 7] (full!)
    /// assert_eq!(min.val, 3);
    /// assert_eq!(min.pos, 1);
    /// ```
    pub fn push(&mut self, val: V, pos: usize) -> Option<MinElem<V>> {
        // Remove strictly larger elements from the back (they'll never be the minimum)
        // For ties, keep the earlier element (smaller position) - leftmost minimum
        while self.deque.back().is_some_and(|back| back.val > val) {
            self.deque.pop_back();
        }

        // Add the new element
        self.deque.push_back(MinElem { pos, val });

        // Remove elements from the front that fall outside the window
        while self
            .deque
            .front()
            .is_some_and(|front| pos >= front.pos + self.w)
        {
            self.deque.pop_front();
        }

        self.current_pos += 1;

        // Return the minimum if the window is full
        if self.current_pos >= self.w {
            self.deque.front().copied()
        } else {
            None
        }
    }

    /// Reset the tracker to its initial state.
    ///
    /// This is useful when processing multiple sequences without reallocating.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::operations::sliding_min::SlidingMin;
    ///
    /// let mut tracker = SlidingMin::new(3).unwrap();
    /// tracker.push(5, 0);
    /// tracker.push(3, 1);
    ///
    /// tracker.reset();
    /// assert_eq!(tracker.push(7, 0), None); // Window not full after reset
    /// ```
    pub fn reset(&mut self) {
        self.deque.clear();
        self.current_pos = 0;
    }

    /// Get the current window size.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::operations::sliding_min::SlidingMin;
    ///
    /// let tracker = SlidingMin::<u64>::new(10).unwrap();
    /// assert_eq!(tracker.window_size(), 10);
    /// ```
    pub fn window_size(&self) -> usize {
        self.w
    }
}

/// Create a sliding window minimum iterator from any iterator.
///
/// This is a convenience function that wraps an iterator and produces
/// the sliding window minima.
///
/// # Arguments
///
/// * `iter` - An iterator over values
/// * `w` - Window size
///
/// # Returns
///
/// An iterator over `MinElem<V>` representing the minimum in each window.
///
/// # Example
///
/// ```
/// use biometal::operations::sliding_min::sliding_min_iter;
///
/// let values = vec![3, 1, 4, 1, 5, 9, 2, 6];
/// let minima: Vec<_> = sliding_min_iter(values.into_iter(), 3)
///     .unwrap()
///     .collect();
///
/// // Windows: [3,1,4], [1,4,1], [4,1,5], [1,5,9], [5,9,2], [9,2,6]
/// assert_eq!(minima[0].val, 1); // min([3,1,4]) = 1
/// assert_eq!(minima[1].val, 1); // min([1,4,1]) = 1
/// assert_eq!(minima[2].val, 1); // min([4,1,5]) = 1
/// ```
pub fn sliding_min_iter<V, I>(
    iter: I,
    w: usize,
) -> Result<impl Iterator<Item = MinElem<V>>>
where
    V: Ord + Copy,
    I: Iterator<Item = V>,
{
    let mut tracker = SlidingMin::new(w)?;

    Ok(iter
        .enumerate()
        .filter_map(move |(pos, val)| tracker.push(val, pos)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_sliding_min() {
        let mut tracker = SlidingMin::new(3).unwrap();

        // Window: [5] (not full)
        assert_eq!(tracker.push(5, 0), None);

        // Window: [5, 3] (not full)
        assert_eq!(tracker.push(3, 1), None);

        // Window: [5, 3, 7] (full!)
        let min1 = tracker.push(7, 2).unwrap();
        assert_eq!(min1.val, 3);
        assert_eq!(min1.pos, 1);

        // Window: [3, 7, 2]
        let min2 = tracker.push(2, 3).unwrap();
        assert_eq!(min2.val, 2);
        assert_eq!(min2.pos, 3);

        // Window: [7, 2, 8]
        let min3 = tracker.push(8, 4).unwrap();
        assert_eq!(min3.val, 2);
        assert_eq!(min3.pos, 3);
    }

    #[test]
    fn test_window_size_one() {
        let mut tracker = SlidingMin::new(1).unwrap();

        let min1 = tracker.push(5, 0).unwrap();
        assert_eq!(min1.val, 5);
        assert_eq!(min1.pos, 0);

        let min2 = tracker.push(3, 1).unwrap();
        assert_eq!(min2.val, 3);
        assert_eq!(min2.pos, 1);
    }

    #[test]
    fn test_invalid_window_size() {
        assert!(SlidingMin::<u64>::new(0).is_err());
    }

    #[test]
    fn test_monotonic_decreasing() {
        // Test case where values are monotonically decreasing
        let mut tracker = SlidingMin::new(3).unwrap();

        assert_eq!(tracker.push(9, 0), None);
        assert_eq!(tracker.push(7, 1), None);

        let min1 = tracker.push(5, 2).unwrap();
        assert_eq!(min1.val, 5);

        let min2 = tracker.push(3, 3).unwrap();
        assert_eq!(min2.val, 3);
    }

    #[test]
    fn test_monotonic_increasing() {
        // Test case where values are monotonically increasing
        let mut tracker = SlidingMin::new(3).unwrap();

        assert_eq!(tracker.push(1, 0), None);
        assert_eq!(tracker.push(3, 1), None);

        let min1 = tracker.push(5, 2).unwrap();
        assert_eq!(min1.val, 1);
        assert_eq!(min1.pos, 0);

        let min2 = tracker.push(7, 3).unwrap();
        assert_eq!(min2.val, 3);
        assert_eq!(min2.pos, 1);
    }

    #[test]
    fn test_all_equal() {
        // Test case where all values are equal (tie-breaking by position)
        let mut tracker = SlidingMin::new(3).unwrap();

        assert_eq!(tracker.push(5, 0), None);
        assert_eq!(tracker.push(5, 1), None);

        let min1 = tracker.push(5, 2).unwrap();
        assert_eq!(min1.val, 5);
        assert_eq!(min1.pos, 0); // Leftmost position wins

        let min2 = tracker.push(5, 3).unwrap();
        assert_eq!(min2.val, 5);
        assert_eq!(min2.pos, 1);
    }

    #[test]
    fn test_reset() {
        let mut tracker = SlidingMin::new(3).unwrap();

        assert_eq!(tracker.push(5, 0), None);
        assert_eq!(tracker.push(3, 1), None);
        tracker.push(7, 2).unwrap();

        tracker.reset();

        assert_eq!(tracker.push(9, 0), None);
        assert_eq!(tracker.push(1, 1), None);

        let min = tracker.push(4, 2).unwrap();
        assert_eq!(min.val, 1);
        assert_eq!(min.pos, 1);
    }

    #[test]
    fn test_sliding_min_iter() {
        let values = vec![3, 1, 4, 1, 5, 9, 2, 6];
        let minima: Vec<_> = sliding_min_iter(values.into_iter(), 3)
            .unwrap()
            .collect();

        assert_eq!(minima.len(), 6); // 8 - 3 + 1 = 6 windows

        // Windows: [3,1,4], [1,4,1], [4,1,5], [1,5,9], [5,9,2], [9,2,6]
        assert_eq!(minima[0].val, 1); // min([3,1,4]) = 1
        assert_eq!(minima[1].val, 1); // min([1,4,1]) = 1
        assert_eq!(minima[2].val, 1); // min([4,1,5]) = 1
        assert_eq!(minima[3].val, 1); // min([1,5,9]) = 1
        assert_eq!(minima[4].val, 2); // min([5,9,2]) = 2
        assert_eq!(minima[5].val, 2); // min([9,2,6]) = 2
    }

    #[test]
    fn test_window_size_getter() {
        let tracker = SlidingMin::<u64>::new(10).unwrap();
        assert_eq!(tracker.window_size(), 10);
    }

    // Property-based test: Compare against naive O(w) implementation
    #[cfg(test)]
    mod properties {
        use super::*;
        use proptest::prelude::*;

        fn naive_sliding_min<V: Ord + Copy>(values: &[V], w: usize) -> Vec<MinElem<V>> {
            if w == 0 || values.len() < w {
                return Vec::new();
            }

            (0..=values.len() - w)
                .map(|i| {
                    let window = &values[i..i + w];
                    let (min_idx, &min_val) = window
                        .iter()
                        .enumerate()
                        .min_by_key(|&(_, &v)| v)
                        .unwrap();
                    MinElem {
                        val: min_val,
                        pos: i + min_idx,
                    }
                })
                .collect()
        }

        proptest! {
            #[test]
            fn prop_matches_naive_implementation(
                values in prop::collection::vec(0u32..100, 1..100),
                w in 1usize..20,
            ) {
                let w = w.min(values.len()); // Ensure w <= values.len()

                let naive_result = naive_sliding_min(&values, w);

                let optimized_result: Vec<_> = sliding_min_iter(values.iter().copied(), w)
                    .unwrap()
                    .collect();

                prop_assert_eq!(naive_result, optimized_result);
            }

            #[test]
            fn prop_minimum_is_in_window(
                values in prop::collection::vec(0u32..100, 1..100),
                w in 1usize..20,
            ) {
                let w = w.min(values.len());

                let minima: Vec<_> = sliding_min_iter(values.iter().copied(), w)
                    .unwrap()
                    .collect();

                for (i, min_elem) in minima.iter().enumerate() {
                    // Check that the minimum position is within the window
                    prop_assert!(min_elem.pos >= i);
                    prop_assert!(min_elem.pos < i + w);

                    // Check that the value matches the position
                    prop_assert_eq!(values[min_elem.pos], min_elem.val);

                    // Check that it's actually the minimum
                    let window = &values[i..i + w];
                    let actual_min = window.iter().copied().min().unwrap();
                    prop_assert_eq!(min_elem.val, actual_min);
                }
            }
        }
    }
}
