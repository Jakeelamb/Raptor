#![allow(
    clippy::cmp_owned,
    clippy::collapsible_if,
    clippy::doc_overindented_list_items,
    clippy::extra_unused_lifetimes,
    clippy::for_kv_map,
    clippy::format_in_format_args,
    clippy::get_first,
    clippy::implicit_saturating_sub,
    clippy::io_other_error,
    clippy::let_and_return,
    clippy::lines_filter_map_ok,
    clippy::manual_div_ceil,
    clippy::manual_hash_one,
    clippy::manual_is_multiple_of,
    clippy::manual_pattern_char_comparison,
    clippy::manual_range_contains,
    clippy::manual_repeat_n,
    clippy::manual_strip,
    clippy::module_inception,
    clippy::needless_borrows_for_generic_args,
    clippy::needless_lifetimes,
    clippy::needless_range_loop,
    clippy::needless_update,
    clippy::new_without_default,
    clippy::ptr_arg,
    clippy::question_mark,
    clippy::redundant_closure,
    clippy::should_implement_trait,
    clippy::sliced_string_as_bytes,
    clippy::too_many_arguments,
    clippy::trim_split_whitespace,
    clippy::type_complexity,
    clippy::unnecessary_cast,
    clippy::unnecessary_literal_unwrap,
    clippy::unwrap_or_default,
    clippy::useless_vec,
    clippy::write_literal
)]

pub mod accel;
pub mod cli;
pub mod dist;
pub mod eval;
pub mod gpu;
pub mod graph;
pub mod io;
pub mod kmer;
pub mod pipeline;
pub mod polish;
pub mod quant;
pub mod stats;
pub mod visualize;
