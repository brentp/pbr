set -eu

target=x86_64-unknown-linux-musl
cargo clean 
cross build --target=$target --release
ls -lhd ./target/$target/release/pbr
