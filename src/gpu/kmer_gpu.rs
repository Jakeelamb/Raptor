use ocl::{ProQue, Buffer, flags};
use std::fs;

pub struct GpuKmerCounter {
    pro_que: ProQue,
    max_reads: usize,
}

impl GpuKmerCounter {
    pub fn new(_k: usize, max_reads: usize) -> Self {
        let kernel_src = fs::read_to_string("src/gpu/kernels/count_kmers.cl").expect("Kernel missing");
        let pro_que = ProQue::builder()
            .src(kernel_src)
            .dims(max_reads)
            .build().unwrap();

        GpuKmerCounter { pro_que, max_reads }
    }

    pub fn count(&self, reads: &Vec<String>, k: usize) -> Vec<u32> {
        // Flatten reads
        let sequence_data: Vec<u8> = reads.iter().flat_map(|s| s.bytes()).collect();

        // Offset positions
        let mut offsets = Vec::with_capacity(reads.len() + 1);
        let mut pos = 0;
        offsets.push(pos as i32);
        for r in reads {
            pos += r.len();
            offsets.push(pos as i32);
        }

        let sequences_buf = Buffer::<u8>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_ONLY)
            .len(sequence_data.len())
            .copy_host_slice(&sequence_data)
            .build().unwrap();

        let offsets_buf = Buffer::<i32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_ONLY)
            .len(offsets.len())
            .copy_host_slice(&offsets)
            .build().unwrap();

        let counts_buf = Buffer::<u32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_WRITE_ONLY)
            .len(65536)
            .build().unwrap();

        let kernel = self.pro_que.kernel_builder("count_kmers")
            .arg(&sequences_buf)
            .arg(&offsets_buf)
            .arg(&counts_buf)
            .arg(k as u32)
            .arg(reads.len() as u32)
            .build().unwrap();

        unsafe { kernel.enq().unwrap(); }

        let mut result = vec![0u32; 65536];
        counts_buf.read(&mut result).enq().unwrap();

        result
    }
}
