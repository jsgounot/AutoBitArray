from aba import AutoBitArray

try :
    from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

except ImportError :
    pass

class Record() :

    def __init__(self, identifier, sequence, name="", description="", char_size=3) :
        self.identifier = identifier
        self.name = name
        self.description = description

        if isinstance(sequence, AutoBitArray) :
            self.sequence = sequence
        else :
            self.sequence = AutoBitArray(sequence, char_size)

    def get_name(self) :
        return self.name if self.name else self.identifier

    def __getitem__(self, key) :
        sub_mba = self.sequence[key]
        return Record(self.identifier, sub_mba, self.name, self.description)

    def __str__(self) :
        return "Record:%s" %self.identifier

    def __repr__(self) :
        return "Record:%s" %self.identifier

    def __len__(self) :
        return len(self.sequence)

    def to_string(self) :
        return "".join(self.sequence.decode())


class MultipleSeqAlignement() :

    def __init__(self, records) :
        self.records = []
        for record in records :
            if not isinstance(record, Record) :
                raise ValueError("records must be composed of Record object only")
            self.records.append(record)

        if not self.records :
            raise ValueError("Empty records")

        if len(set(len(record.sequence) for record in self.records)) != 1 :
            raise ValueError("Record should have the same length")

    def __len__(self) :
        return len(self.records[0].sequence)

    def __getitem__(self, key) :
        return MultipleSeqAlignement([record.__getitem__(key)
        for record in self.records])

    def to_biopython(self, code=unambiguous_dna) :
        return MultipleSeqAlignment(SeqRecord(Seq(record.to_string(), code),
        id=record.identifier, name = record.name, description=record.description)
        for record in self.records)

    def save(self, outfile, kind, mode="w", **kwargs) :
        save_method = {"fasta" : self.to_fasta,
        "nexus" : self.to_nexus,
        "phylip" : self.to_phylip,
        "structure" : self.to_structure}

        fun = save_method.get(kind, None)
        if fun is None :
            lm = list(save_method.keys())
            raise ValueError("Format not found in the list %s" %lm)

        fun(outfile, mode, **kwargs)

    def get_names_same_size(self, maxsize=None, addspace=False, exception_same=True, space_add=0) :
        larger = max(len(record.get_name()) for record in self.records)
        names = [record.get_name() + " " * (larger - len(record.get_name()) + space_add)
        for record in self.records]

        if maxsize is not None :
            newnames = [record[:maxsize] for record in names]
            if len(set(newnames)) != len(names) and exception_same :
                raise ValueError("Duplicated names after truncation")
            names = newnames

        if addspace and maxsize :
            names = [name + " " * (maxsize - len(name)) for name in names]

        return names

    @staticmethod
    def safe_names(names) :
        WHITESPACE = '\t\n'
        PUNCTUATION = '()[]{}\,;:=*\'"`+-<>'
        return ["'%s'" %name if set(name).intersection(set(WHITESPACE + PUNCTUATION)) else name
        for name in names]

    def to_fasta(self, outfile, mode, blocksize=60) :

        with open(outfile, mode) as f :
            for record in self.records :
                f.write(">%s\n" %record.get_name())
                for i in range(0, len(record), blocksize) :
                    f.write("%s\n" %record[i:i+blocksize].to_string())

    def to_nexus(self, outfile, mode, dtype="dna", missing="?", gap="-", interleaves=True,
                 blocksize=70, safename=True) :

        with open(outfile, mode) as f :
            nchar = len(self)
            ntax = len(self.records)

            f.write("#NEXUS\n")
            f.write("begin data;\n")
            f.write("\t" + "dimensions ntax=%i nchar=%i;\n" %(ntax, nchar))
            f.write("\t" + "format datatype=%s missing=%s gap=%s interleave;\n" %(dtype, missing, gap))
            f.write("matrix\n")

            names = self.get_names_same_size()
            names = MultipleSeqAlignement.safe_names(names) if safename else names

            if interleaves :
                for i in range(0, len(self), blocksize) :
                    for idx, name in enumerate(names) :
                        f.write("%s %s\n" %(name, self.records[idx][i:i+blocksize].to_string()))
                    f.write("\n")

            else :
                for idx, name in enumerate(names) :
                    f.write("%s %s\n" %(name, self.records[idx].to_string()))

            f.write(";\nend;")

    def to_phylip(self, outfile, mode, strict=True, blocksize=10, nblockline=5) :

        names = self.get_names_same_size(10, addspace=True) if strict else self.get_names_same_size()
        max_len_name = max(len(name) for name in names)
        first_loop = True
        total_block_size = blocksize * nblockline

        with open(outfile, mode) as f :
            nchar = len(self)
            ntax = len(self.records)

            f.write(" %i %i\n" %(ntax, nchar))
            for i in range(0, len(self), total_block_size) :
                for idx, name in enumerate(names) :
                    record = self.records[idx][i:i+total_block_size]

                    if first_loop :
                        f.write(name)

                    else :
                        f.write(" " * max_len_name)

                    for j in range(nblockline) :
                        subindex = blocksize * j
                        subrec = record[subindex:subindex+blocksize].to_string()
                        f.write(" %s" %subrec)

                    f.write("\n")

                if first_loop :
                    first_loop = False

                f.write("\n")

    def to_structure(self, outfile, mode, positions, code=None, fun_name=None) :

        code = code or {"A" : "1", "T" : "2", "G" : "3", "C" : "4"}

        positions = list(positions)
        positions = [str(position - positions[idx - 1]) if idx != 0 else "-1"
                     for idx, position in enumerate(positions)]

        nunique = set()
        with open(outfile, mode) as f :
            f.write("\t".join(positions) + "\n")

            for record in sorted(self.records, key = lambda x : x.get_name()) :
                name = record.get_name()
                if fun_name : name = fun_name(name)
                nunique.add(name)

                encoded_seq = (code.get(base, "-9") for base in record.to_string())
                f.write(str(name) + "\t" + "\t".join(encoded_seq) + "\n")                

        print ("Structure writer : %s" %outfile)
        print ("NUMINDS : %i" %len(nunique))
        print ("NUMLOCI : %i" %len(self))
        print ("MISSING : -9")