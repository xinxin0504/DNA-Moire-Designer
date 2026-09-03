# Moiré sequence-I/O extension scope

Base application: official `cadnano2 2.4.13`.

The Windows companion changes only legacy JSON I/O and positional-file launch:

1. read and restore top-level `scaffold_sequences` records at their 5′ anchor;
2. save all non-empty scaffold sequences through ordinary Save/Save As;
3. preserve unknown top-level Moiré Designer metadata fields;
4. preserve/regenerate scaffold colours, `lattice`, and `num_bases`;
5. open a JSON passed as the first positional command-line argument.

No Designer UI, routing, SST generation, staple generation, capture generation,
sequence design, or analysis code is added to the companion cadnano program.
