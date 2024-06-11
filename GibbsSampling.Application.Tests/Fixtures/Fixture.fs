module Fixture

open BioFSharp

open CompositeVector.Functions
open PositionMatrix.Functions

type Fixture() =
    
    static member  Pseudocount = 1.
    
    static member DNABases =
        [|
            Nucleotides.Nucleotide.A
            Nucleotides.Nucleotide.T
            Nucleotides.Nucleotide.G
            Nucleotides.Nucleotide.C
            Nucleotides.Nucleotide.Gap
        |]

    static member AminoAcids =
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

    static member ShortProfileSequences = "CACGTG"
    
    static member LongProfileSequences =
        [|
            "GTGGCTGCACCACGTGTATGC"
            "ACATCGCATCACGTGACCAGT"
            "CCTCGCACGTGGTGGTACAGT"
            "CTCGTTAGGACCATCACGTGA"
        |]

    static member ComplexProfileSequences =
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

    static member FrequencyCompositeVector =        
        BioArray.ofNucleotideString Fixture.ShortProfileSequences
        |> createFCVOf
       
    static member ProbabilityCompositeVector =        
        createPCVOf Fixture.FrequencyCompositeVector

    static member PositionProbabilityMatrix =        
        BioArray.ofNucleotideString Fixture.ShortProfileSequences
        |> createPFMOf
        |> createPPMOf
        |> normalizePPM Fixture.ShortProfileSequences.Length Fixture.DNABases Fixture.Pseudocount

    static member PositionWeightMatrix =
        new PositionWeightMatrix(Array2D.length2 Fixture.PositionProbabilityMatrix.Matrix)
