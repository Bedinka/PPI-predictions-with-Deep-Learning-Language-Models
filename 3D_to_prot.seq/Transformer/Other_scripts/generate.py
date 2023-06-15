def make_small_data(infile, outfile):
    f = open(infile)
    fo = open(outfile, "w")
    i = 0
    for line in f:
        print(i)
        if len(line) > 500:
            start = (len(line) - 500) // 2
            end = start + 500
            fo.write(line[start:end] + "\n")
        else:
            fo.write(line)
        i += 1
        if i > 150000:
            break
    f.close()
    fo.close()

if __name__ == "__main__":
    make_small_data("./data/train_fasta.90.txt", "./data/train_fasta.train12.txt")
    make_small_data("./data/train_structure.90.txt", "./data/train_structure.train12.txt")