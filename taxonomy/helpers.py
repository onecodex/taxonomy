import hashlib


def md5_for_file(filename, block_size=2**20):
    md5 = hashlib.md5()
    with open(filename, mode='rb') as f:
        while True:
            data = f.read(block_size)
            if not data:
                break
            md5.update(data)
    return md5.hexdigest()
