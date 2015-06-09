# How to

Generate table.

    $ make gen
    $ ./tablegen TABLESIZE > output.dat

TABLESIZE would be 6~10(256B ~ 1KB) for Epiphany, since Epiphany only have 8KB/bank local mem.

Include table data for exp() code.
