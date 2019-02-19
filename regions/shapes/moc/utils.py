def trailing_zeros(x):
    bits = 0
    # convention for x == 0 => return 0
    if x == 0:
        return 0

    if not (x & 0xFFFFFFFF):
        bits += 32
        x >>= 32
    if not (x & 0xFFFF):
        bits += 16
        x >>= 16
    if not (x & 0xFF):
        bits += 8
        x >>= 8
    if not (x & 0xF):
        bits += 4
        x >>= 4
    if not (x & 0x3):
        bits += 2
        x >>= 2
    if not (x & 1):
        bits += 1
        x >>= 1

    return bits
