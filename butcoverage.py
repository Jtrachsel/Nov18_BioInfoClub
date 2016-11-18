__author__ = 'julian'
from sys import argv
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.IUPACData import ambiguous_dna_values
import pandas as pd
import re
from Bio.Phylo.Applications import FastTreeCommandline
from statistics import mode
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Align.Applications import MuscleCommandline



def length_filter(sequences, minlength=None, maxlength=None):
    lengths = []
    for seq in sequences:
        lengths.append(len(seq.seq))

    if minlength is None:
        minlength = min(lengths)
    if maxlength is None:
        maxlength = max(lengths)
    goodseqs_min = []
    goodseqs_max = []
    for seq in sequences:
        if len(seq.seq) <= maxlength:
            goodseqs_max.append(seq)
        else:
            # print("sequence {} is longer than the maximum length and has been removed".format(seq.id))
            continue
    for seq in goodseqs_max:
        if len(seq.seq) >= minlength:
            goodseqs_min.append(seq)
        else:
            # print("sequence {} is shorter than the minimum length and has been removed".format(seq.id))
            continue
    return goodseqs_min


def codon_align(protein_alignment, dna_sequences):

    codon_positions = []

    for seq in protein_alignment:
        codons = []
        codon_count = 0
        for position in seq.seq:
            if position != "-":
                codons.append(codon_count)
                codon_count += 3
            else:
                codon_count += 3
        codon_positions.append(codons)

    DNA_sequence_index = 0
    Codon_list_index = 0

    DNA_align_length = len(protein_alignment[1].seq) * 3
    Codon_DNA_alignment = []

    for x in codon_positions:
        blank_alignment = "-" * DNA_align_length
        blank_alignment = MutableSeq(blank_alignment)
        DNA_align = blank_alignment
        DNA_position_index = 0
        for y in x:
            DNA_align[y:y+3] = dna_sequences[DNA_sequence_index].seq[DNA_position_index:(DNA_position_index+3)]
            DNA_position_index += 3

        DNA_align = DNA_align.toseq()
        DNA_align = SeqRecord(DNA_align)
        DNA_align.id = dna_sequences[DNA_sequence_index].id
        DNA_align.name = dna_sequences[DNA_sequence_index].name
        DNA_align.description = dna_sequences[DNA_sequence_index].description
        DNA_align.dbxrefs = dna_sequences[DNA_sequence_index].dbxrefs
        DNA_sequence_index += 1
        Codon_DNA_alignment.append(DNA_align)
        Codon_list_index += 1
    return Codon_DNA_alignment


def unique_seqs(sequences):
    """returns a list of SeqRecord objects with redundant sequences removed"""
    unique_records = []
    checksum_container = []
    for seq in sequences:
        checksum = seguid(seq.seq)
        if checksum not in checksum_container:
            checksum_container.append(checksum)
            unique_records.append(seq)
    return unique_records


def unique_ids(sequences):
    """returns a list of SeqRecord objects with redundant ids renamed"""
    unique_records = []
    checksum_container = []
    redundant_id_count = 0
    for seq in sequences:
        checksum = seguid(seq.id)
        if checksum not in checksum_container:
            checksum_container.append(checksum)
            unique_records.append(seq)
        else:
            print("repeated id detected, adding '.{}' suffix".format(redundant_id_count))
            seq.id = "{}.{}".format(seq.id, redundant_id_count)
            unique_records.append(seq)
            redundant_id_count += 1
    return unique_records


def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""

    return SeqRecord(seq=nuc_record.seq.translate(cds=False), id=nuc_record.id, description="translation using def table")


def alignment_screen(protein_alignment, protein_seqs, dna_seqs, cutoff):
    aln = protein_alignment
    alignment_length = (len(aln[0]))
    alignment_iterator = range(0, alignment_length)

    iterator_count = 0
    list_of_percent_gaps = []
    for x in alignment_iterator:
        column = aln[:, iterator_count]
        number_of_gaps = column.count('-')
        percent_gaps = number_of_gaps/len(column)
        list_of_percent_gaps.append(percent_gaps)
        iterator_count += 1
    # print(list_of_percent_gaps)

    counter = 0
    good_columns = []
    bad_columns = []
    for entry in list_of_percent_gaps:
        if entry < cutoff:
            good_columns.append(counter)
            counter += 1
        else:
            bad_columns.append(counter)
            counter += 1

    prot_DNA_dict = {}
    for sequence in protein_alignment:
        position_counter = 0
        residue_counter = 0
        prot_DNA_dict[sequence.id] = {}
        for position in sequence.seq:
            if position != "-":
                prot_DNA_dict[sequence.id][position_counter] = residue_counter
                residue_counter += 1
            position_counter += 1

    residues_removed = {}
    for entry in prot_DNA_dict:
        # print(entry)
        residues_removed[entry] = []
        for position in prot_DNA_dict[entry]:
            if position in bad_columns:
                residues_removed[entry].append(prot_DNA_dict[entry][position])

    #print(residues_removed)
    nucleotides_removed = {}
    for entry in residues_removed:
        nucleotides_removed[entry] = []
        for y in residues_removed[entry]:
            y = [y*3, (y*3) + 1, (y*3) + 2]
            nucleotides_removed[entry].append(y)
        nucleotides_removed[entry] = sum(nucleotides_removed[entry], [])  # flattens the lists of lists into just lists
    #for entry in nucleotides_removed:
        #print("sequence ID: {}\nNucleotides removed: {}".format(entry, nucleotides_removed[entry]))

    prot_dict = {}
    dna_dict = {}

    for sequence in protein_seqs:
        prot_dict[sequence.id] = sequence.seq

    for sequence in dna_seqs:
        dna_dict[sequence.id] = sequence.seq

    edited_prots = []

    for seq in protein_seqs:
        blank_seq = []
        position_counter = 0
        # print(len(seq.seq))
        #print("id = {}".format(seq.id))
        for position in seq.seq:
            if position_counter in residues_removed[seq.id]:

                #print("RESIDUE #{} REMOVED".format(position_counter))
                position_counter += 1
            else:
                position_counter += 1
                blank_seq.append(position)
        new_seq = "".join(blank_seq)
        seq.seq = Seq(new_seq)
        edited_prots.append(seq)

    edited_DNAs = []

    for seq in dna_seqs:
        blank_seq = []
        position_counter = 0
        for position in seq.seq:
            if position_counter in nucleotides_removed[seq.id]:
                position_counter += 1
            else:
                position_counter += 1
                blank_seq.append(position)
        new_seq = "".join(blank_seq)
        seq.seq = Seq(new_seq)
        edited_DNAs.append(seq)

    return edited_prots, edited_DNAs


def find_hit_regions(primer, alignment): #this one is for all the sequences in the alignment
    '''this is currently super inefficient... It basically does the work of primer_coverage() for every single possible
    frame in a sliding window for every sequence... If I'm ok with this I should just have this function return the
    number of mismatches for the positions which best match...  If I do that then I could have the amplicon length be
    something that was returned as well.....hmmm very tempting... I think I should do this.  what else besides amplicon
    length would this allow me to do?  I could also have it output potential mispriming sites, and then the amplicon
     length for the misprimed sites.... I could include a condition where it would print a warning if mispriming
     is likely, output a spreadsheet that tells you what sequences are likely to misprime, how big the amplicon
     for the mispriming would be...  But this mispriming would only be for these particular sequences that you are
     tyring to amplify, A much more liekly source of mispriming would just be other random genomic DNA.  A metagenome
     might be a good thing to run this, but that would really take a long time.....'''

    alignment_len = len(alignment[0])
    primer_length = len(primer)
    number_of_frames = (alignment_len - primer_length) + 1
    range_of_frames = range(0, number_of_frames)
    list_of_indexes = []
    first_indexes = []
    last_indexes = []
    frame_indexes = {}
    for frame in range_of_frames:
        frame_indexes[frame] = {}
        frame_indexes[frame]["first"] = frame
        frame_indexes[frame]["last"] = frame + primer_length

    hit_regions = {}
    for seq in alignment:
        sequences = {}
        for frame in frame_indexes:
            sequence = seq[frame_indexes[frame]["first"]:frame_indexes[frame]["last"]]
            #print(sequence)
            sequences[frame] = sequence

        number_mismatches = {}
        for key in sequences:
            number_mismatches[key] = 0
            for count, position in enumerate(sequences[key].upper()):
                #print(count, position)
                if position not in ambiguous_dna_values[primer[count]]:
                    number_mismatches[key] += 1
        indexes = frame_indexes[min(number_mismatches, key=number_mismatches.get)]
        hit_regions[seq.id] = indexes
        #print("number of sequences checked: {}".format(len(hit_regions)))
        #print("Percent complete: {}".format(len(hit_regions)/len(alignment)))
    #hit_regions = set(hit_regions)
    #print(hit_regions)

    starting = []
    ending = []
    for key in hit_regions:
        #print(key)
        starting.append(hit_regions[key]["first"])
        ending.append(hit_regions[key]["last"])
    #print(starting)
    #print(ending)
    starting = mode(starting)
    ending = mode(ending)
    return starting, ending


def primer_coverage(FWDprimer, REVprimer, FWDregion, REVregion):
    """ Returns a pandas dataframe with the number of mismatches for each primer in their corresponding region.
    It also attempts to pull metadata from the sequences (organism, and gene definition) this works if the sequences
    are from fungene.  Row names are sequence IDs. """
    number_mismatches_FWD = {}
    melt_temp_FWD = {}

    for seq in FWDregion:
        number_mismatches_FWD[seq.id] = 0
        melt_temp_FWD[seq.id] = mt.Tm_NN(seq.seq)
        for count, position in enumerate(seq.seq.upper()):
            if position not in ambiguous_dna_values[FWDprimer[count]]:
                number_mismatches_FWD[seq.id] += 1

    REVprimer = REVprimer.reverse_complement()

    number_mismatches_REV = {}
    melt_temp_REV = {}
    for seq in REVregion:
        number_mismatches_REV[seq.id] = 0
        melt_temp_REV[seq.id] = mt.Tm_NN(seq.seq)
        for count, position in enumerate(seq.seq.upper()):
            if position not in ambiguous_dna_values[REVprimer[count]]:
                number_mismatches_REV[seq.id] += 1

    mismatches_fwd_rev = {}      # there is probably a better way to name this dict....
    for key in number_mismatches_FWD:
        mismatches_fwd_rev[key] = {}
        mismatches_fwd_rev[key]["FWD_mismatch"] = number_mismatches_FWD[key]
        mismatches_fwd_rev[key]["REV_mismatch"] = number_mismatches_REV[key]
        mismatches_fwd_rev[key]["FWD_Melt_temp"] = melt_temp_FWD[key]
        mismatches_fwd_rev[key]["REV_Melt_temp"] = melt_temp_REV[key]

    mismatches_df = pd.DataFrame.from_dict(mismatches_fwd_rev, 'index')

    # making the metadata whatnot

    metadata = {}
    org_regex = re.compile(r"organism=.*,")
    def_regex = re.compile(r"definition=.*,")
    fun_regex = re.compile(r"function=.*")
    for seq in FWDregion:
        metadata[seq.id] = {}
        organism = re.search(org_regex, seq.description).group(0)
        organism = (organism[9:(len(organism)-1)])  #gets rid of the "organism=" and "," that my regex pulled out
        definition = re.search(def_regex, seq.description).group(0)
        definition = (definition[11:(len(definition)-1)])  #gets rid of the "definition=" that my regex pulled out
        function = re.search(fun_regex, seq.description).group(0)
        function = function[9:(len(function))]
        metadata[seq.id]["organism"] = organism
        metadata[seq.id]["definition"] = definition
        metadata[seq.id]["function"] = function

    metadata_df = pd.DataFrame.from_dict(metadata, 'index')

    final_df = pd.merge(mismatches_df, metadata_df, left_index=True, right_index=True)
    return final_df


########################  Parameters  #####################
gap_cutoff = 0.95
''' if any position in the protein alignment has a gap in greater than this percent of all sequences, this position will
 be removed from the protein alignment, and the corresponding codon will be removed from the DNA codon alignment'''
###########################################################

dna = SeqIO.parse(argv[1], "fasta") # reads in the unaligned DNA sequences downloaded from FunGene
dna = list(dna) # turns the generator returned from SeqIO.parse into a list we can iterate through more than once
print("you input {} DNA sequences".format(len(dna)))
print("filtering sequences based on length, maxlength = 1600, minlength = 1200...")
dna = length_filter(dna, 1200, 1700) # removes all sequences less than 1200bp or more than 1700bp
print("There are {} DNA sequences after length filter".format(len(dna)))
dna = unique_seqs(dna)  # removes redundant identical sequences
print("there are {} unique DNA sequences".format(len(dna)))
dna = unique_ids(dna) #



prot = []
for seq in dna:
    seq = make_protein_record(seq)
    prot.append(seq)


unique_dna_handle = "uniqueDNA.{}".format(argv[1])
protein_handle = "protein.{}".format(argv[1])
protein_align_handle = "align.{}".format(protein_handle)
unique_prot_dna_handle = "unique.prot.DNA.{}".format(argv[1])


unique_proteins = unique_seqs(prot)     # this little block removes redundant protein sequences
unique_protein_ids = []
for seq in unique_proteins:
    unique_protein_ids.append(seq.id)


unique_prots_dna = []
# generates the final DNA seqs list based on the the protein seqs that made the cut
for seq in dna:
    if seq.id in unique_protein_ids:
        unique_prots_dna.append(seq)


SeqIO.write(unique_proteins, protein_handle, "fasta")
print("protein sequences output to {}".format(protein_handle))

SeqIO.write(unique_prots_dna, unique_prot_dna_handle, "fasta")
print("final DNA sequences output to {}".format(unique_prot_dna_handle))

print("there are {} sequences in the analysis".format(len(unique_prots_dna)))
print("aligning sequences with Clustal Omega...")


#cline = ClustalOmegaCommandline(infile=protein_handle, outfile=protein_align_handle, auto=True)
#cline()


cline = MuscleCommandline(input=protein_handle, out=protein_align_handle, maxiters=2)
cline()


#THESE ARE NECESSARY BECAUSE MUSCLE SUCKS AND REARRANGES YOUR SEQUENCES AND THEIR OPTION TO STOP THAT BEHAVIOR IS BROKEN
align = SeqIO.parse(protein_align_handle, "fasta")
sortedList = [f for f in sorted(align, key=lambda x : x.id)]
SeqIO.write(sortedList, protein_align_handle, 'fasta')


dna = SeqIO.parse(unique_prot_dna_handle, 'fasta')
sortedList2 = [f for f in sorted(dna, key=lambda x : x.id)]
SeqIO.write(sortedList2, unique_prot_dna_handle, 'fasta')
unique_prots_dna = SeqIO.parse(unique_prot_dna_handle, 'fasta')
unique_prots_dna = list(unique_prots_dna)


print("Done! Alignment output to {}".format(protein_align_handle))

prot_align = AlignIO.read(protein_align_handle, "fasta")


dna_codon_align = codon_align(prot_align, unique_prots_dna)
dna_codon_align_handle = "codon.align.{}".format(argv[1])
print("generating DNA codon alignment...")
SeqIO.write(dna_codon_align, dna_codon_align_handle, "fasta")
print("Done! DNA codon alignment output to {}".format(dna_codon_align_handle))

print("screening protein alignment, removing residues in columns that have >95% gaps")
screened_prot, screened_dna = alignment_screen(prot_align, unique_proteins, unique_prots_dna, gap_cutoff)
screened_dna_handle = "screened.DNA.{}".format(argv[1])
screened_prot_handle = "screened.prots.{}".format(argv[1])
SeqIO.write(screened_prot, screened_prot_handle, "fasta")

SeqIO.write(screened_dna, screened_dna_handle, "fasta")

print("realigning screened protein residues...")
screened_protein_aligned_handle = "screened.prot.align.{}".format(argv[1])
#cline1 = ClustalOmegaCommandline(infile=screened_prot_handle, outfile=screened_protein_aligned_handle, auto=True)
#cline1()

cline1 = MuscleCommandline(input=screened_prot_handle, out=screened_protein_aligned_handle, maxiters=2)
cline1()


#THESE ARE NECESSARY BECAUSE MUSCLE SUCKS AND REARRANGES YOUR SEQUENCES AND THEIR OPTION TO STOP THAT BEHAVIOR IS BROKEN
align = SeqIO.parse(screened_protein_aligned_handle, "fasta")
sortedList = [f for f in sorted(align, key=lambda x : x.id)]
SeqIO.write(sortedList, screened_protein_aligned_handle, 'fasta')
screened_prot_align = AlignIO.read(screened_protein_aligned_handle, "fasta")


dna = SeqIO.parse(screened_dna_handle, 'fasta')
sortedList2 = [f for f in sorted(dna, key=lambda x : x.id)]
SeqIO.write(sortedList2, screened_dna_handle, 'fasta')
screened_dna = SeqIO.parse(screened_dna_handle, 'fasta')
screened_dna = list(screened_dna)
print("Done! Screened protein alignment output to {}".format(screened_protein_aligned_handle))

print("regenerating DNA codon alignment with screened sequences")

screened_codon_align = codon_align(screened_prot_align, screened_dna)
screened_codon_align_handle = "screened.codon.align.{}".format(argv[1])
SeqIO.write(screened_codon_align, screened_codon_align_handle, "fasta")
print("DONE!!!! Screened DNA codon alignment output to {}".format(screened_codon_align_handle))

screened_codon_align = MultipleSeqAlignment(screened_codon_align) # do I need this??

vitalF1 = Seq("CAGCTNGGYATYGGNGS", IUPAC.ambiguous_dna)
vitalF2 = Seq("GGWATWGGMGSYATGCC", IUPAC.ambiguous_dna)
vitalF3 = Seq("GHATYGGNGSTATGCC", IUPAC.ambiguous_dna)

vitalR1 = Seq("AARTCCANYTGNCCVCC", IUPAC.ambiguous_dna)
vitalR2 = Seq("AARTCCANYTGNCCVCC", IUPAC.ambiguous_dna)
vitalR3 = Seq("AAGTCWAAYTGWCCRCC", IUPAC.ambiguous_dna)
vitalR1rc = vitalR1.reverse_complement()
vitalR2rc = vitalR2.reverse_complement()
vitalR3rc = vitalR3.reverse_complement()

BCoATscrF = Seq("GCNGANCATTTCACNTGGAAYWSNTGGCAYATG", IUPAC.ambiguous_dna)
BCoATscrR = Seq("CCTGCCTTTGCAATRTCNACRAANGC", IUPAC.ambiguous_dna)
BCoATscrRrc = BCoATscrR.reverse_complement()

MybutFWD = Seq("CARYTNGGNATYGGNGGNATSCC", IUPAC.ambiguous_dna)
MybutREV = Seq("TGTCCGCCNGYNCCRSWRAT", IUPAC.ambiguous_dna)
MybutREVrc = MybutREV.reverse_complement()


print("Finding optimal binding regions for each primer.  This will take a very long time because my code is lazy")


BCoATscrF_indexes = find_hit_regions(BCoATscrF, screened_codon_align)
print('one primer done, out of 10 total')
BCoATscrR_indexes = find_hit_regions(BCoATscrRrc, screened_codon_align)
print('two...done with the flint primers')
FlintFWDregion = (screened_codon_align[:, BCoATscrF_indexes[0]:BCoATscrF_indexes[1]])
#print(FlintFWDregion)
FlintREVregion = (screened_codon_align[:, BCoATscrR_indexes[0]:BCoATscrR_indexes[1]])
#print(FlintREVregion)
vitalF1_indexes = find_hit_regions(vitalF1, screened_codon_align)
print('three...')
vitalF2_indexes = find_hit_regions(vitalF2, screened_codon_align)
print('four')
vitalF3_indexes = find_hit_regions(vitalF3, screened_codon_align)
print('five')
vitalR1_indexes = find_hit_regions(vitalR1rc, screened_codon_align)
print('six')
#vitalR2_indexes = find_hit_regions(vitalR2rc, screened_codon_align)
print('seven')
#vitalR3_indexes = find_hit_regions(vitalR3rc, screened_codon_align)
print('eight....done with the Vital primers')
vitalF1region = (screened_codon_align[:, vitalF1_indexes[0]:vitalF1_indexes[1]])
vitalF2region = (screened_codon_align[:, vitalF2_indexes[0]:vitalF2_indexes[1]])
vitalF3region = (screened_codon_align[:, vitalF3_indexes[0]:vitalF3_indexes[1]])
vitalR1region = (screened_codon_align[:, vitalR1_indexes[0]:vitalR1_indexes[1]])
vitalR2region = (screened_codon_align[:, vitalR1_indexes[0]:vitalR1_indexes[1]])
vitalR3region = (screened_codon_align[:, vitalR1_indexes[0]:vitalR1_indexes[1]])


#print(vitalF1region)
#print(vitalF2region)
#print(vitalF3region)

#print(vitalR1region)
#print(vitalR2region)
#print(vitalR3region)


MyFWDindexes = find_hit_regions(MybutFWD, screened_codon_align)
print('nine')
MyREVindexes = find_hit_regions(MybutREVrc, screened_codon_align)
print('ten!')
MyFWDregion = (screened_codon_align[:, MyFWDindexes[0]:MyFWDindexes[1]])
MyREVregion = (screened_codon_align[:, MyREVindexes[0]:MyREVindexes[1]])
#print(MyFWDregion)
#print(MyREVregion)


print('running primer_coverage()')
vitalF1_R1_df = primer_coverage(vitalF1, vitalR1, vitalF1region, vitalR1region)
vitalF1_R1_df = vitalF1_R1_df[["FWD_mismatch", "REV_mismatch", "FWD_Melt_temp", "REV_Melt_temp", "organism", "definition", "function"]]
vitalF1_R1_df.columns = ["vitalF1_MM", "vitalR1_MM", "vitalF1_Melt", "vitalR1_Melt", "organism", "definition", "function"]
vitalF1_R1_df.to_csv('vitalF1_R1_df.primermismatch.csv', sep="\t")
vitalF1_R1_df = vitalF1_R1_df.drop(["definition", "function", "organism"], 1)

vitalF2_R2_df = primer_coverage(vitalF2, vitalR2, vitalF2region, vitalR2region)
vitalF2_R2_df = vitalF2_R2_df[["FWD_mismatch", "REV_mismatch", "FWD_Melt_temp", "REV_Melt_temp", "organism", "definition", "function"]]
vitalF2_R2_df.columns = ["vitalF2_MM", "vitalR2_MM", "vitalF2_Melt", "vitalR2_Melt", "organism", "definition", "function"]
vitalF2_R2_df.to_csv('vitalF2_R2_df.primermismatch.csv', sep="\t")
vitalF2_R2_df = vitalF2_R2_df.drop(['definition', 'function', 'organism'], 1)

vitalF3_R3_df = primer_coverage(vitalF3, vitalR3, vitalF3region, vitalR3region)
vitalF3_R3_df = vitalF3_R3_df[["FWD_mismatch", "REV_mismatch", "FWD_Melt_temp", "REV_Melt_temp", "organism", "definition", "function"]]
vitalF3_R3_df.columns = ["vitalF3_MM", "vitalR3_MM", "vitalF3_Melt", "vitalR3_Melt", "organism", "definition", "function"]
vitalF3_R3_df.to_csv('vitalF3_R3_df.primermismatch.csv', sep="\t")
vitalF3_R3_df = vitalF3_R3_df.drop(['definition', 'function', 'organism'], 1)

Myprimers_df = primer_coverage(MybutFWD, MybutREV, MyFWDregion, MyREVregion)
Myprimers_df = Myprimers_df[["FWD_mismatch", "REV_mismatch", "FWD_Melt_temp", "REV_Melt_temp", "organism", "definition", "function"]]
Myprimers_df.columns = ["but672FWD_MM", "but1031REV_MM", "but672FWD_Melt", "but1031REV_Melt", "organism", "definition", "function"]
Myprimers_df.to_csv('Myprimers_df.primermismatch.csv', sep="\t")
Myprimers_df = Myprimers_df.drop(['definition', 'function', 'organism'], 1)


Flint_df = primer_coverage(BCoATscrF, BCoATscrR, FlintFWDregion, FlintREVregion)
Flint_df = Flint_df[["FWD_mismatch", "REV_mismatch", "FWD_Melt_temp", "REV_Melt_temp", "organism", "definition", "function"]]
Flint_df.columns = ["BCoATscrF_MM", "BCoATscrR_MM", "BCoATscrF_Melt", "BCoATscrR_Melt", "organism", "definition", "function"]
Flint_df.to_csv('Flint_df.primermismatch.csv', sep="\t")

combined_df = pd.merge(vitalF1_R1_df, vitalF2_R2_df, left_index=True, right_index=True, how='outer')
combined_df = pd.merge(combined_df, vitalF3_R3_df, left_index=True, right_index=True, how='outer')


combined_df = pd.merge(combined_df, Myprimers_df, left_index=True, right_index=True, how='outer')
combined_df = pd.merge(combined_df, Flint_df, left_index=True, right_index=True, how='outer')

combined_df.to_csv("finalbutcoverage.csv", sep="\t")


print("building tree.....")

fasttree_cline = FastTreeCommandline("FastTreeMP", input=protein_align_handle, out="Tree1.nwk")

print('using FastTree with the command: {}'.format(fasttree_cline))
fasttree_cline()

print('done!  Thanks!!')



