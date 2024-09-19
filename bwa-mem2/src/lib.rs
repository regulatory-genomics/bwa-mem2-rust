use std::{ffi::{CStr, CString}, path::Path};
use noodles::sam::{self, header::record::value::map::ReferenceSequence};
use noodles::sam::header::record::value::Map;
use noodles::fastq;
use bstr::ByteSlice;

#[derive(Debug)]
pub struct IndexError(String);

/// A BWA reference object to perform alignments to.
/// Must be loaded from a BWA index created with `bwa index`
#[derive(Debug)]
pub struct FMIndex {
    fm_index: bwa_mem2_sys::FMI_search,
    ref_bytes: Vec<u8>,
    contig_names: Vec<String>,
    contig_lengths: Vec<usize>,
}

impl FMIndex {
    pub fn new<P1: AsRef<Path>, P2: AsRef<Path>>(fasta: P1, location: P2) -> Result<Self, IndexError> {
        let fasta_file = CString::new(fasta.as_ref().to_str().unwrap()).unwrap();
        let prefix = CString::new(location.as_ref().to_str().unwrap()).unwrap();
        unsafe {
            bwa_mem2_sys::bwa_idx_build(fasta_file.as_ptr(), prefix.as_ptr());
        }
        Self::read(location)
    }

    /// Load a BWA reference from disk. Pass the fasta filename of the
    /// original reference as `path`
    pub fn read<P: AsRef<Path>>(path: P) -> Result<Self, IndexError> {
        let idx_file = CString::new(path.as_ref().to_str().unwrap()).unwrap();
        let mut contig_names = Vec::new();
        let mut contig_lengths = Vec::new();
        let ref_file = path.as_ref().to_string_lossy() + ".0123";
        let ref_bytes = std::fs::read(ref_file.as_ref())
            .expect(&format!("Could not read reference file {}", ref_file));

        unsafe {
            let mut idx = bwa_mem2_sys::FMI_search::new(idx_file.as_ptr());
            idx.load_index();

            let num_contigs = (*((*(idx._base.idx)).bns)).n_seqs;
            for i in 0..num_contigs as isize {
                let name = CStr::from_ptr((*((*((*(idx._base.idx)).bns)).anns.offset(i))).name);
                let sz = (*((*((*(idx._base.idx)).bns)).anns.offset(i))).len;

                let name_string = name.to_owned().into_string().unwrap();
                contig_names.push(name_string);
                contig_lengths.push(sz as usize)
            }
            Ok(Self { fm_index: idx, contig_lengths, contig_names, ref_bytes })
        }

    }

    pub fn create_sam_header(&self) -> sam::Header {
        let ref_seqs = self.contig_names.iter().zip(self.contig_lengths.iter()).map(|(name, len)|
            (bstr::BString::from(name.as_str()), Map::<ReferenceSequence>::new(std::num::NonZeroUsize::try_from(*len).unwrap()))
        ).collect();
        sam::Header::builder().set_reference_sequences(ref_seqs).build()
    }

    pub fn filename(&self) -> &str {
        unsafe {
            CStr::from_ptr(&self.fm_index.file_name as *const i8).to_str().unwrap()
        }
    }
}

impl Drop for FMIndex {
    fn drop(&mut self) {
        unsafe {
            self.fm_index.destruct();
        }
    }
}

/// BWA opts object. Currently only default opts are enabled
pub struct AlignerOpts {
    opts: bwa_mem2_sys::mem_opt_t,
}

impl Default for AlignerOpts {
    fn default() -> Self {
        Self::new()
    }
}

impl AlignerOpts {
    /// Create a `Bwaopts` object with default BWA parameters
    pub fn new() -> Self {
        let ptr = unsafe { bwa_mem2_sys::mem_opt_init() };
        let opts = unsafe { *ptr };
        unsafe { libc::free(ptr as *mut libc::c_void) };
        Self { opts }
    }

    pub fn get_actual_chunk_size(&self) -> usize {
        (self.opts.chunk_size * self.opts.n_threads as i64).try_into().unwrap()
    }

    pub fn set_n_threads(mut self, n_threads: usize) -> Self {
        self.opts.n_threads = n_threads as i32;
        self
    }

    /// Set alignment scores
    pub fn set_scores(
        mut self,
        matchp: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        self.opts.a = matchp;
        self.opts.b = mismatch;
        self.opts.o_del = gap_open;
        self.opts.o_ins = gap_open;
        self.opts.e_del = gap_extend;
        self.opts.e_ins = gap_extend;

        unsafe {
            bwa_mem2_sys::bwa_fill_scmat(matchp, mismatch, self.opts.mat.as_mut_ptr());
        }
        self
    }

    /// Set clipping score penalties
    pub fn set_clip_scores(mut self, clip5: i32, clip3: i32) -> Self {
        self.opts.pen_clip5 = clip5;
        self.opts.pen_clip3 = clip3;
        self
    }

    /// Set unpaired read penalty
    pub fn set_unpaired(mut self, unpaired: i32) -> Self {
        self.opts.pen_unpaired = unpaired;
        self
    }

    /// Mark shorter splits as secondary
    pub fn set_no_multi(mut self) -> Self {
        self.opts.flag |= 0x10; // MEM_F_NO_MULTI
        self
    }
}

/// Paired-end statistics structure used by BWA to score paired-end reads
pub struct PairedEndStats {
    inner: [bwa_mem2_sys::mem_pestat_t; 4],
}

impl Default for PairedEndStats {
    fn default() -> Self {
        Self::simple(200.0, 100.0, 35, 600)
    }
}

impl PairedEndStats {
    /// Generate a 'simple' paired-end read structure that standard forward-reverse
    /// pairs as created by TruSeq, Nextera, or Chromium Genome sample preparations.
    pub fn simple(avg: f64, std: f64, low: i32, high: i32) -> PairedEndStats {
        let pe_stat_null = || bwa_mem2_sys::mem_pestat_t {
            failed: 1,
            low: 0,
            high: 0,
            avg: 0.0,
            std: 100.0,
        };

        let pes = [
            pe_stat_null(),
            bwa_mem2_sys::mem_pestat_t {
                failed: 0,
                low,
                high,
                avg,
                std,
            },
            pe_stat_null(),
            pe_stat_null(),
        ];

        PairedEndStats { inner: pes }
    }
}

/// A BWA aligner. Carries everything required to align
/// reads to a reference and generate BAM records.
pub struct BurrowsWheelerAligner {
    index: FMIndex,
    opts: AlignerOpts,
    header: sam::Header,
    pe_stats: PairedEndStats,
}

struct WorkerWrapper {
    ptr: *mut bwa_mem2_sys::worker_t,
    n_threads: i32,
}

impl WorkerWrapper {
    fn new(num_reads: i32, n_threads: i32) -> Self {
        let ptr = unsafe { bwa_mem2_sys::new_worker_t(num_reads, n_threads) };
        Self { ptr, n_threads }
    }
}

impl Drop for WorkerWrapper {
    fn drop(&mut self) {
        unsafe {
            bwa_mem2_sys::destroy_worker_t(self.ptr, self.n_threads);
        }
    }
}

impl BurrowsWheelerAligner {
    pub fn new(
        index: FMIndex,
        opts: AlignerOpts,
        pe_stats: PairedEndStats,
    ) -> Self {
        let header = index.create_sam_header();
        Self { index, opts, header, pe_stats }
    }
    pub fn get_sam_header(&self) -> sam::Header {
        self.header.clone()
    }

    pub fn align_reads(&mut self, records: &mut [fastq::Record]) -> impl ExactSizeIterator<Item = sam::Record> {
        //unsafe { bwa_mem2_sys::bwa_verbose = 0; }
        let mut seqs: Vec<_> = records.iter_mut().enumerate().map(|(i, fq)| new_bseq1_t(i, fq)).collect();
        let mut opts = self.opts.opts.clone();
        let num_reads = seqs.len().try_into().unwrap();
        unsafe {
            let worker = WorkerWrapper::new(num_reads, self.opts.opts.n_threads);
            (*worker.ptr).ref_string = self.index.ref_bytes.as_ptr() as *mut u8;
            (*worker.ptr).nthreads = self.opts.opts.n_threads.try_into().unwrap();
            (*worker.ptr).nreads = num_reads;
            (*worker.ptr).fmi = &mut self.index.fm_index as *mut bwa_mem2_sys::FMI_search;

            //let index = self.index.fm_index;
            bwa_mem2_sys::mem_process_seqs(
                &mut opts,
                0,
                num_reads,
                seqs.as_mut_ptr(), 
                self.pe_stats.inner.as_ptr(),
                worker.ptr,
            );
            seqs.into_iter().map(|seq| {
                let sam = CStr::from_ptr(seq.sam).to_str().unwrap().as_bytes().try_into().unwrap();
                libc::free(seq.sam as *mut libc::c_void);
                sam
            })
        }
    }
}

fn new_bseq1_t(id: usize, fq: &mut fastq::Record) -> bwa_mem2_sys::bseq1_t {
    bwa_mem2_sys::bseq1_t {
        name: CString::new(fq.name().as_bytes()).unwrap().into_raw(),
        id: id.try_into().unwrap(),
        l_seq: fq.sequence().len() as i32,
        seq: fq.sequence_mut().as_mut_ptr() as *mut i8,
        qual: fq.quality_scores_mut().as_mut_ptr() as *mut i8,
        comment: std::ptr::null_mut(),
        sam: std::ptr::null_mut(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{fs, io::BufReader};
    use noodles::sam;
    use noodles::fastq;
    use flate2;
    use itertools::Itertools;

    #[test]
    fn test_index() {
        let fasta = "tests/data/ecoli.fa.gz";
        let location = "tests/temp";
        let _ = fs::remove_dir_all(location);
        let idx = FMIndex::new(fasta, location).unwrap();
        println!("{:?}", idx);
    }

    #[test]
    fn test_align() {
        let location = "tests/temp";
        let idx = FMIndex::read(location).unwrap();
        let opt = AlignerOpts::default();
        let mut aligner = BurrowsWheelerAligner::new(idx, opt, PairedEndStats::default());
        let header = aligner.get_sam_header();

        let mut writer = sam::io::Writer::new(std::fs::File::create("out.sam").unwrap());
        let reader = flate2::read::MultiGzDecoder::new(std::fs::File::open("tests/test.fq.gz").unwrap());
        let mut reader = fastq::Reader::new(BufReader::new(reader));
        reader.records().map(|x| x.unwrap()).chunks(100000).into_iter().for_each(|chunk| {
            let mut records: Vec<_> = chunk.collect();
            for sam in aligner.align_reads(&mut records) {
                writer.write_record(&header, &sam).unwrap();
            }
        });
    }
}