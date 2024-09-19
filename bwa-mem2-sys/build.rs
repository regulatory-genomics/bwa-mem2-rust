use std::env;
use std::fs;
use std::path::Path;

fn build_libsafestring() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let obj_dir = Path::new(&out_dir).join("obj");

    // Create the output directories
    fs::create_dir_all(&obj_dir).unwrap();

    // Define source and header files
    let src_dir = "ext/bwa-mem2/ext/safestringlib/safeclib";
    let include_dir = "ext/bwa-mem2/ext/safestringlib/include";
    let c_files = [
        "abort_handler_s.c", "stpcpy_s.c", "strlastsame_s.c", "ignore_handler_s.c", 
        "stpncpy_s.c", "strljustify_s.c", "memcmp16_s.c", "strcasecmp_s.c", 
        "strncat_s.c", "memcmp32_s.c", "strcasestr_s.c", "strncpy_s.c", 
        "memcmp_s.c", "strcat_s.c", "strnlen_s.c", "memcpy16_s.c", 
        "strcmpfld_s.c", "strnterminate_s.c", "memcpy32_s.c", "strcmp_s.c", 
        "strpbrk_s.c", "memcpy_s.c", "strcpyfldin_s.c", "strprefix_s.c", 
        "memmove16_s.c", "strcpyfldout_s.c", "strremovews_s.c", "memmove32_s.c", 
        "strcpyfld_s.c", "strspn_s.c", "memmove_s.c", "strcpy_s.c", 
        "strstr_s.c", "mem_primitives_lib.c", "strcspn_s.c", "strtok_s.c", 
        "strfirstchar_s.c", "strtolowercase_s.c", "memset16_s.c", "strfirstdiff_s.c", 
        "strtouppercase_s.c", "memset32_s.c", "strfirstsame_s.c", "strzero_s.c", 
        "memset_s.c", "strisalphanumeric_s.c", "wcpcpy_s.c", "memzero16_s.c", 
        "strisascii_s.c", "wcscat_s.c", "memzero32_s.c", "strisdigit_s.c", 
        "wcscpy_s.c", "memzero_s.c", "strishex_s.c", "wcsncat_s.c", 
        "strislowercase_s.c", "wcsncpy_s.c", "safe_mem_constraint.c", 
        "strismixedcase_s.c", "wcsnlen_s.c", "strispassword_s.c", "wmemcmp_s.c", 
        "safe_str_constraint.c", "strisuppercase_s.c", "wmemcpy_s.c", 
        "strlastchar_s.c", "wmemmove_s.c", "snprintf_support.c", "strlastdiff_s.c", 
        "wmemset_s.c"
    ];

    let mut build = cc::Build::new();

    // Set compiler flags
    build
        .warnings(false)  // Disable warnings
        .include(include_dir)
        .flag("-fstack-protector-strong")
        .flag("-fPIE")
        .flag("-fPIC")
        .flag("-O2")
        .define("_FORTIFY_SOURCE", "2")
        .flag("-Wformat")
        .flag("-Wformat-security")
        .flag("-z noexecstack")
        .flag("-z relro")
        .flag("-z now");

    // Compile each C file and add it to the build
    for file in &c_files {
        build.file(Path::new(src_dir).join(file));
    }

    // Compile and generate the static library
    build.compile("libsafestring.a");
}

fn main() {
    build_libsafestring();

    // Set up the build process for libbwa.a
    let mut build = cc::Build::new();

    // Determine the architecture flags based on environment variables
    let arch = env::var("CARGO_CFG_TARGET_FEATURE").unwrap_or_default();
    println!("cargo:warning=arch: {}", arch);
    let arch_flags = if arch.contains("sse4.1") {
        vec!["-msse", "-msse2", "-msse3", "-mssse3", "-msse4.1"]
    } else if arch.contains("sse4.2") {
        vec!["-msse", "-msse2", "-msse3", "-mssse3", "-msse4.1", "-msse4.2"]
    } else if arch.contains("avx2") {
        vec!["-mavx2"]
    } else if arch.contains("avx512") {
        vec!["-mavx512bw"]
    } else if arch.contains("native") {
        vec!["-march=native"]
    } else {
        vec!["-msse", "-msse2", "-msse3", "-mssse3", "-msse4.1"]
    };

    for f in arch_flags.into_iter() {
        build.flag(f);
    }

    // Set compiler flags
    build
        .warnings(false)  // Disable warnings
        .cpp(true)
        .flag("-g")
        .flag("-O3")
        .flag("-fpermissive")
        .flag("-fPIC")
        .define("ENABLE_PREFETCH", None)
        .define("V17", Some("1"))
        .define("MATE_SORT", Some("0"))
        .include("ext/bwa-mem2/src")
        .include("ext/bwa-mem2/ext/safestringlib/include");

    // List of source files
    let src_dir = "ext/bwa-mem2/src/";
    let files = [
        "fastmap.cpp", "bwtindex.cpp", "utils.cpp", "memcpy_bwamem.cpp",
        "kthread.cpp", "kstring.cpp", "ksw.cpp", "bntseq.cpp", "bwamem.cpp",
        "profiling.cpp", "bandedSWA.cpp", "FMI_search.cpp", "read_index_ele.cpp",
        "bwamem_pair.cpp", "kswv.cpp", "bwa.cpp", "bwamem_extra.cpp", "kopen.cpp",
    ];

    // Add files to build
    build.file("wrapper.cpp");
    for file in &files {
        build.file(Path::new(src_dir).join(file));
    }

    // Compile to static library
    build.compile("libbwa.a");

    println!("cargo:rustc-link-search=native={}", env::var("OUT_DIR").unwrap());

    // Link the pre-built safestringlib
    //println!("cargo:rustc-link-lib=static=safestring");

    //println!("cargo:rustc-link-lib=static=bwa");

    // Link the zlib library
    println!("cargo:rustc-link-lib=z");

    // Link the C++ standard library
    //println!("cargo:rustc-link-lib=stdc++");

}
