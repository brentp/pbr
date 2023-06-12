set -eu

#target=x86_64-unknown-linux-gnu
#export RUSTFLAGS="-C target-feature=-crt-static -C relocation-model=pic"
#cargo test --release --target $target && cargo build --release --target $target
#ls -lh ./target/$target/release/pbr
#exit

target=x86_64-unknown-linux-musl
cargo clean 
cross build --target=$target --release
ls -lhd ./target/$target/release/pbr
