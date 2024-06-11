module Fixture

open BioFSharp

type Fixture() =

    let longProfileSequences =
        [|
            "GTGGCTGCACCACGTGTATGC"
            "ACATCGCATCACGTGACCAGT"
            "CCTCGCACGTGGTGGTACAGT"
            "CTCGTTAGGACCATCACGTGA"
        |]

    let shortProfileSequences =
        [|
            "CACGTG"
            "CACGTG"
            "CACGTG"
            "CACGTG"
        |]

    let thyminFocusedProfileSequences =
        [|
            "GTGGCTGCACCACGTGTATGCCACGTG"
            "ACATCGCATCACGTGACCAGTTAGTTG"
            "CCTCGCACGTGGTGGTACAGTCGTACG"
            "GCATAAAGGACCATCACGTGAAGCTGC"
            "TTTTTTTTTTTTTTTTTTTTTTTTTTT"
        |]

    let veryLongProfileSequences =
        [|
            "GTAAGTACAGAAAGCCACAGAGTACCATCTAGGAAATTAACATTATACTAACTTTCTACATCGTTGATACTTATGCGTATACATTCATATA"
            "AGACAGAGTCTAAAGATTGCATTACAAGAAAAAAGTTCTCATTACTAACAAGCAAAATGTTTTGTTTCTCCTTTTA"
            "GTATGTTCATGTCTCATTCTCCTTTTCGGCTCCGTTTAGGTGATAAACGTACTATATTGTGAAAGATTATTTACTAACGACACATTGAAG*"
            "GCATGTGTGCTGCCCAAGTTGAGAAGAGATACTAACAAAATGACCGCGGCTCTCAAAAATAATTGACGAGCTTACGGTGATACGCTTACCG"
            "GTATGTTTGACGAGAATTGCTAGTGTGCGGGAAACTTTGCTACCTTTTTTGGTGCGATGCAACAGGTTACTAATATGTAATACTTCAG"
            "TTTCAAGATTAACCACATCTGCTAACTTTCTCCCTATGCTTTTACTAACAAAATTATTCTCACTCCCCGATATTGA"
            "GTAAGTATCCAGATTTTACTTCATATATTTGCCTTTTTCTGTGCTCCGACTTACTAACATTGTATTCTCCCCTTCTTCATTTTAG"
            "GTATGCATAGGCAATAACTTCGGCCTCATACTCAAAGAACACGTTTACTAACATAACTTATTTACATAG"
            "GTATGTAGTAGGGAAATATATCAAAGGAACAAAATGAAAGCTATGTGATTCCGTAATTTACGAAGGCAAATTACTAACATTGAAATACGGG"
            "GTATGTTACTATTTGGAGTTTCATGAGGCTTTTCCCGCCGTAGATCGAACCCAATCTTACTAACAGAGAAAGGGCTTTTTCCCGACCATCA"
            "TATGTAATGATATATTATGAAGTAAGTTCCCCAAAGCCAATTAACTAACCGAATTTTAATCTGCACTCATCATTAG"
            "GTATGTTCATAATGATTTACATCGGAATTCCCTTTGATACAAGAAAACTAACGGGTATCGTACATCAATTTTTGAAAAAAGTCAAGTACTA"
            "GTATGTATATTTTTGACTTTTTGAGTCTCAACTACCGAAGAGAAATAAACTACTAACGTACTTTAATATTTATAG"
            "TTTCGACGCGAATAGACTTTTTCCTTCTTACAGAACGATAATAACTAACATGACTTTAACAG"
        |]

    let dnaBases =
        [|
            Nucleotides.Nucleotide.A
            Nucleotides.Nucleotide.T
            Nucleotides.Nucleotide.G
            Nucleotides.Nucleotide.C
            Nucleotides.Nucleotide.Gap
        |]

    let aminoAcids =
        [|
            AminoAcidSymbols.AminoAcidSymbol.Ala
            AminoAcidSymbols.AminoAcidSymbol.Arg
            AminoAcidSymbols.AminoAcidSymbol.Asn
            AminoAcidSymbols.AminoAcidSymbol.Asp
            AminoAcidSymbols.AminoAcidSymbol.Asx
            AminoAcidSymbols.AminoAcidSymbol.Cys
            AminoAcidSymbols.AminoAcidSymbol.Xle
            AminoAcidSymbols.AminoAcidSymbol.Gln
            AminoAcidSymbols.AminoAcidSymbol.Glu
            AminoAcidSymbols.AminoAcidSymbol.Glx
            AminoAcidSymbols.AminoAcidSymbol.Gly
            AminoAcidSymbols.AminoAcidSymbol.His
            AminoAcidSymbols.AminoAcidSymbol.Ile
            AminoAcidSymbols.AminoAcidSymbol.Leu
            AminoAcidSymbols.AminoAcidSymbol.Lys
            AminoAcidSymbols.AminoAcidSymbol.Met
            AminoAcidSymbols.AminoAcidSymbol.Phe
            AminoAcidSymbols.AminoAcidSymbol.Pro
            AminoAcidSymbols.AminoAcidSymbol.Pyl
            AminoAcidSymbols.AminoAcidSymbol.Sel
            AminoAcidSymbols.AminoAcidSymbol.Ser
            AminoAcidSymbols.AminoAcidSymbol.Thr
            AminoAcidSymbols.AminoAcidSymbol.Trp
            AminoAcidSymbols.AminoAcidSymbol.Val
        |]
