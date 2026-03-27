use std::convert::TryFrom;
use std::io::{self, ErrorKind, Read, Write};
use std::path::Path;

use tempfile::{Builder, NamedTempFile};

const BRANCH_READ_SPOOL_MAGIC: [u8; 4] = *b"BRSP";
const BRANCH_READ_SPOOL_VERSION: u32 = 1;

#[derive(Debug)]
pub struct BranchReadSpool {
    file: NamedTempFile,
}

impl BranchReadSpool {
    pub fn create(temp_dir: Option<&Path>) -> io::Result<Self> {
        let file = match temp_dir {
            Some(dir) => Builder::new()
                .prefix("raptor-branch-read-spool")
                .tempfile_in(dir)?,
            None => Builder::new()
                .prefix("raptor-branch-read-spool")
                .tempfile()?,
        };

        let mut spool = Self { file };
        spool.write_header()?;
        Ok(spool)
    }

    pub fn append(&mut self, sequence: &[u8]) -> io::Result<()> {
        let len = u32::try_from(sequence.len())
            .map_err(|_| io::Error::new(ErrorKind::InvalidInput, "sequence exceeds u32::MAX"))?;

        self.file.as_file_mut().write_all(&len.to_le_bytes())?;
        self.file.as_file_mut().write_all(sequence)?;
        Ok(())
    }

    pub fn replay_sequences<F>(&self, mut on_sequence: F) -> io::Result<()>
    where
        F: FnMut(&[u8]) -> io::Result<()>,
    {
        let mut reader = self.open_reader()?;
        self.read_header(&mut reader)?;

        loop {
            let mut len_buf = [0u8; 4];
            match reader.read_exact(&mut len_buf) {
                Ok(()) => {}
                Err(err) if err.kind() == ErrorKind::UnexpectedEof => break,
                Err(err) => return Err(err),
            }

            let len = u32::from_le_bytes(len_buf) as usize;
            let mut sequence = vec![0u8; len];
            reader.read_exact(&mut sequence)?;
            on_sequence(&sequence)?;
        }

        Ok(())
    }

    pub fn replay_into_batches<F>(&self, batch_size: usize, mut on_batch: F) -> io::Result<()>
    where
        F: FnMut(&[Vec<u8>]) -> io::Result<()>,
    {
        let batch_size = batch_size.max(1);
        let mut batch = Vec::with_capacity(batch_size);
        self.replay_sequences(|sequence| {
            batch.push(sequence.to_vec());
            if batch.len() >= batch_size {
                on_batch(&batch)?;
                batch.clear();
            }
            Ok(())
        })?;

        if !batch.is_empty() {
            on_batch(&batch)?;
        }

        Ok(())
    }

    fn open_reader(&self) -> io::Result<std::io::BufReader<std::fs::File>> {
        Ok(std::io::BufReader::new(self.file.reopen()?))
    }

    fn write_header(&mut self) -> io::Result<()> {
        let file = self.file.as_file_mut();
        file.write_all(&BRANCH_READ_SPOOL_MAGIC)?;
        file.write_all(&BRANCH_READ_SPOOL_VERSION.to_le_bytes())?;
        Ok(())
    }

    fn read_header<R: Read>(&self, reader: &mut R) -> io::Result<()> {
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if magic != BRANCH_READ_SPOOL_MAGIC {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "branch read spool magic mismatch",
            ));
        }

        let mut version = [0u8; 4];
        reader.read_exact(&mut version)?;
        if u32::from_le_bytes(version) != BRANCH_READ_SPOOL_VERSION {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "branch read spool version mismatch",
            ));
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::BranchReadSpool;

    #[test]
    fn branch_read_spool_round_trips_sequences() {
        let mut spool = BranchReadSpool::create(None).expect("create spool");
        let reads = [b"ACGTACGT".as_slice(), b"TTTT".as_slice(), b"".as_slice()];

        for read in &reads {
            spool.append(read).expect("append read");
        }

        let mut replayed = Vec::new();
        spool
            .replay_sequences(|sequence| {
                replayed.push(sequence.to_vec());
                Ok(())
            })
            .expect("replay sequences");

        assert_eq!(
            replayed,
            reads.iter().map(|read| read.to_vec()).collect::<Vec<_>>()
        );
    }
}
