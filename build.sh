set -eu

#target=x86_64-unknown-linux-gnu
#export RUSTFLAGS="-C target-feature=-crt-static -C relocation-model=pie"
#cargo test --release --target $target 
#RUSTFLAGS="-C target-feature=-crt-static -C relocation-model=pie" cargo build --release --target $target
#ls -lh ./target/$target/release/pbr
#exit

target=x86_64-unknown-linux-musl
cargo clean 
cross build --target=$target --release
ls -lhd ./target/$target/release/pbr
