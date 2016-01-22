#!/usr/bin/env python
from pbtools.pbtranscript.io import GFF
from collections import defaultdict
from pbtools.pbtranscript.counting import compare_junctions
from pbcore.io.FastqIO import FastqReader, FastqWriter
from csv import DictReader, DictWriter


def can_merge(m, r1, r2, internal_fuzzy_max_dist):
    if m == 'subset':
        r1, r2 = r2, r1 #  rotate so r1 is always the longer one
    if m == 'super' or m == 'subset':
        n2 = len(r2.ref_exons)
        if r1.strand == '+':
            return abs(r1.ref_exons[-1].start - r2.ref_exons[-1].start) <= internal_fuzzy_max_dist and \
                r1.ref_exons[-n2].start <= r2.ref_exons[0].start < r1.ref_exons[-n2].end
        else:
            return abs(r1.ref_exons[0].end - r2.ref_exons[0].end) <= internal_fuzzy_max_dist and \
                    r1.ref_exons[n2-1].start <= r2.ref_exons[-1].end < r1.ref_exons[n2].end

def filter_out_subsets(recs, internal_fuzzy_max_dist):
    # recs must be sorted by start becuz that's the order they are written
    i = 0
    while i < len(recs)-1:
        j = i + 1
        while j < len(recs):
            if recs[j].start > recs[i].end: 
                break
            recs[i].segments = recs[i].ref_exons
            recs[j].segments = recs[j].ref_exons
            m = compare_junctions.compare_junctions(recs[i], recs[j], internal_fuzzy_max_dist)
            if can_merge(m, recs[i], recs[j], internal_fuzzy_max_dist):
                if m == 'super': # pop recs[j] 
                    recs.pop(j)
                else:
                    recs.pop(i)
                    j += 1
            else:
                j += 1
        i += 1


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_prefix", help="Input prefix")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("--fuzzy_junction", type=int, default=5, help="Fuzzy junction max dist (default: 5bp)")

    args = parser.parse_args()

    #group_filename = args.input_prefix + '.group.txt'
    count_filename = args.input_prefix + '.abundance.txt'
    gff_filename = args.input_prefix + '.gff'
    rep_filename = args.input_prefix + '.rep.fq'

    recs = defaultdict(lambda: [])
    reader = GFF.collapseGFFReader(gff_filename)
    for r in reader:
        assert r.seqid.startswith('PB.')
        recs[int(r.seqid.split('.')[1])].append(r)

    good = []
    f = open(args.output_prefix + '.gff', 'w')
    keys = recs.keys()
    keys.sort()
    for k in recs:
        xxx = recs[k]
        filter_out_subsets(xxx, args.fuzzy_junction)
        for r in xxx:
            GFF.write_collapseGFF_format(f, r)
            good.append(r.seqid)
    f.close()

    # read abundance first
    f = open(count_filename)
    count_header = ''
    while True:
        cur_pos = f.tell()
        line = f.readline()
        if not line.startswith('#'):
            f.seek(cur_pos)
            break
        else:
            count_header += line
    d = dict((r['pbid'], r) for r in DictReader(f, delimiter='\t'))
    for k,v in d.iteritems():
        print k,v
    f.close()

    # write output rep.fq
    f = FastqWriter(args.output_prefix + '.rep.fq')
    for r in FastqReader(rep_filename):
        if r.name.split('|')[0] in good:
           f.writeRecord(r)
    f.close()

    # write output to .abundance.txt
    f = open(args.output_prefix + '.abundance.txt', 'w')
    f.write(count_header)
    writer = DictWriter(f, fieldnames=['pbid','count_fl','count_nfl','count_nfl_amb','norm_fl','norm_nfl','norm_nfl_amb'], \
                        delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for k in good:
        r = d[k]
        writer.writerow(r)
    f.close()


if __name__ == "__main__":
    main()
