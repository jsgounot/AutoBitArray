from bitarray import bitarray
from itertools import product
from copy import deepcopy
from math import log, ceil

class BitDict() :

    def __init__(self, char_size=3) :
        self._dbit = {}
        self.char_size = char_size
        self.init_bit_gen(char_size)

    def init_bit_gen(self, char_size) :
        self._retrieved = []
        self._bit_gen = product('01', repeat=char_size)

    @staticmethod
    def char_size_from_nunique(nunique) :
        # Return the char_size which should be used based on n
        # number of unique elements

        nunique = nunique if isinstance(nunique, int) else len(set(nunique))
        return int(ceil(log(nunique, 2)))

    def encode_element(self, element) :
        # Return the bit representation of a single element
        # If element is not already found in the internal dict
        # A new representation is created and stored
        # Raise an exception if no more representations are available

        try :
            return self._dbit[element]

        except KeyError :
            bit_rep = next(self._bit_gen)
            bit_rep = bitarray(''.join(bit_rep))
            return self._dbit.setdefault(element, bit_rep)

    def replace(self, key, new_key) :
        if key not in self._dbit :
            raise ValueError("Element not found in the dictionnary")

        if new_key in self._dbit :
            raise ValueError("New element found in the dictionnary")            

        else :
            self._dbit[new_key] = self._dbit.pop(key)

    def encode(self, elements) :
        # Encode list of values into one single bitarray

        for element in set(elements) :
            try : 
                self.encode_element(element)
            except StopIteration : 
                raise ValueError("The number of unique elements exceed the char_size possibilities")

        b = bitarray()
        b.encode(self._dbit, elements)
        return b

    def decode(self, seq) :
        # Decode the bitarray

        return seq.decode(self._dbit)


class AutoBitArray() :

    def __init__(self, seq, char_size=3, bit_dict=None) :
        self.bdict = bit_dict if bit_dict else BitDict(char_size)
        self.init_seq(seq)

    def init_seq(self, seq) :
        self._seq = bitarray()
        self.extend(seq)

    def subnew(self, * barrays) :
        # Create a new mba with the same bdict based
        # on zero to multiples bit arrays

        new = AutoBitArray("", bit_dict=self.bdict)
        for barray in barrays :
            new._seq.extend(barray)
        return new

    def copy(self) :
        new = AutoBitArray("", bit_dict=self.bdict)
        new._seq = self._seq.copy()
        return new      

    def __add__(self, other) :
        if not isinstance(other, AutoBitArray) :
            self.extend(other)

        elif self.bdict != other.bdict :
            raise ValueError("Only two mbas with the same BitDict can be summed")
        
        else :
            return self.subnew(self._seq, other._seq)

    def __str__(self) :
        if len(self) > 10 :
            return str(self[:5]) + "..." + str(self[-5:])

        else :
            return "".join(self.bdict.decode(self._seq))

    def __len__(self) :
        return int(len(self._seq) / self.bdict.char_size)

    def __setitem__(self, key, items) :
        key = self.transform_index(key)
        self._seq[key] = self.bdict.encode(items)

    def __getitem__(self, key) :
        key = self.transform_index(key)
        return self.subnew(self._seq[key])

    def replace(self, key, new_key) :
        self.bdict.replace(key, new_key)

    def decode(self) :
        return self.bdict.decode(self._seq)

    def extend(self, items) :
        self._seq.extend(self.bdict.encode(items))

    def transform_index(self, key) :
        if isinstance(key, int) :
            return slice(key * self.bdict.char_size, (key + 1) * self.bdict.char_size, None)
        
        else :
            # slice object
            start, end, step = key.start, key.stop, key.step

            if step is not None :
                raise Exception("Slice with step is currently not implemented")

            values = [start, end, step]
            for idx, value in enumerate(values) :
                if value is not None :
                    values[idx] = value * self.bdict.char_size
            return slice(* values)

