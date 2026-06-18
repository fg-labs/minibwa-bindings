use std::env;
use std::path::PathBuf;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = manifest.join("vendor/minibwa");
    let shim = manifest.join("shim");
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());

    // docs.rs containers lack system zlib and cannot compile the vendored C.
    // Run bindgen only (it needs the vendored headers + libclang, which docs.rs
    // provides) and return early — no C compilation, no link-lib directives.
    if env::var("DOCS_RS").is_ok() {
        generate_bindings(&manifest, &vendor, &shim, &out);
        return;
    }

    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let is_x86 = target_arch == "x86_64";
    let is_aarch64 = target_arch == "aarch64";
    let is_debug = env::var("DEBUG").as_deref() == Ok("true");

    if env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos") {
        env::set_var("MACOSX_DEPLOYMENT_TARGET", "11.0");
    }

    // Object set = LOBJS ∪ {index, bseq, libsais, libsais64, kthread}, excluding
    // the ksw2 SIMD TUs (compiled separately below) and CLI-only objects
    // (main, fastmap, map-main). libsais.c (32-bit) is required because
    // libsais64.c calls libsais_int()/libsais(); kthread.c is kept for the
    // future mb_map_batch threading path (harmless, no main).
    let base_srcs = [
        "kommon.c",
        "kalloc.c",
        "bwt.c",
        "l2bit.c",
        "options.c",
        "seed.c",
        "map-algo.c",
        "lchain.c",
        "align.c",
        "pe.c",
        "cs.c",
        "format.c",
        "index.c",
        "bseq.c",
        "libsais.c",
        "libsais64.c",
        "kthread.c",
    ];

    let mut build = cc::Build::new();
    build.cpp(false).std("c99").include(&vendor).include(&shim);
    if !is_debug {
        build.define("NDEBUG", None);
    }
    for s in base_srcs {
        build.file(vendor.join(s));
    }
    build.file(shim.join("minibwa_shim.c"));

    // ksw2 dispatch + reference impl (all arches).
    build.file(vendor.join("ksw2_extd2_dispatch.c"));
    // 128-bit SSE TUs (NEON via s2n-lite.h on arm64).
    for s in ["ksw2_extz2_sse.c", "ksw2_ll_sse.c"] {
        build.file(vendor.join(s));
    }

    // GCC on aarch64 is stricter than clang about NEON intrinsic vector types in
    // the bundled s2n-lite.h SSE->NEON shim (e.g. passing a uint32x4_t where an
    // int32x4_t is expected). Allow same-size vector conversions, matching what
    // clang does implicitly, so the NEON shim compiles under both compilers.
    if is_aarch64 {
        build.flag_if_supported("-flax-vector-conversions");
    }

    if is_x86 {
        build.flag_if_supported("-msse4.2");
        build.flag_if_supported("-mpopcnt");
        // Wide AVX2 / AVX-512 dual-gap TUs.
        for (tier, flags) in [
            ("ksw_extd2_avx2", &["-mavx2"][..]),
            ("ksw_extd2_avx512", &["-mavx512bw"][..]),
        ] {
            let mut k = cc::Build::new();
            k.cpp(false).std("c99").include(&vendor).include(&shim);
            if !is_debug {
                k.define("NDEBUG", None);
            }
            k.file(vendor.join("ksw2_extd2_wide.c"));
            k.define("ksw_extd2_sse", tier);
            for f in flags {
                k.flag_if_supported(f);
            }
            k.compile(&format!("minibwa-{tier}"));
        }
    }

    if cfg!(feature = "gpl-bwtgen") {
        build.file(vendor.join("QSufSort.c"));
        build.file(vendor.join("bwtgen.c"));
        build.define("USE_GPL", None);
    }
    if cfg!(feature = "openmp") {
        build.flag_if_supported("-fopenmp");
        build.define("LIBSAIS_OPENMP", None);
        println!("cargo:rustc-link-lib=gomp");
    }

    build.warnings(false);
    build.compile("minibwa");

    // ksw2_extd2_sse.c compiled as the reference impl — linked AFTER the main
    // minibwa lib so that extd2_ref_impl (defined here) resolves the forward
    // declaration in ksw2_extd2_dispatch.c which lives in the main lib.
    {
        let mut k = cc::Build::new();
        k.cpp(false).std("c99").include(&vendor).include(&shim);
        if !is_debug {
            k.define("NDEBUG", None);
        }
        k.file(vendor.join("ksw2_extd2_sse.c"));
        k.define("ksw_extd2_sse", "extd2_ref_impl");
        if is_x86 {
            k.flag_if_supported("-msse4.2");
            k.flag_if_supported("-mpopcnt");
        }
        if is_aarch64 {
            k.flag_if_supported("-flax-vector-conversions");
        }
        k.warnings(false);
        k.compile("minibwa-ksw-ref");
    }

    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rerun-if-changed=shim/minibwa_shim.c");
    println!("cargo:rerun-if-changed=shim/minibwa_shim.h");
    println!("cargo:rerun-if-changed=vendor/COMMIT");

    generate_bindings(&manifest, &vendor, &shim, &out);
}

fn generate_bindings(
    manifest: &std::path::Path,
    vendor: &std::path::Path,
    shim: &std::path::Path,
    out: &std::path::Path,
) {
    let _ = manifest; // unused outside docs.rs path but kept for symmetry
    let bindings = bindgen::Builder::default()
        .header(shim.join("minibwa_shim.h").to_str().unwrap())
        .header(vendor.join("minibwa.h").to_str().unwrap())
        .clang_arg(format!("-I{}", vendor.display()))
        .clang_arg(format!("-I{}", shim.display()))
        .allowlist_function("mb_.*")
        .allowlist_type("mb_.*")
        .allowlist_var("MB_.*")
        .layout_tests(true)
        .generate()
        .expect("bindgen failed");
    bindings.write_to_file(out.join("bindings.rs")).unwrap();
}
