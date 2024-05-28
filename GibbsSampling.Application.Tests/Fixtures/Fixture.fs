module Fixture

open BioFSharp

open CompositeVector.Functions
open PositionMatrix.Functions

type Fixture() =
    
    member this.Pseudocount with get () = 1.
    
    member this.DNABases with get () =
        [|
            Nucleotides.Nucleotide.A
            Nucleotides.Nucleotide.T
            Nucleotides.Nucleotide.G
            Nucleotides.Nucleotide.C
            Nucleotides.Nucleotide.Gap
        |]

    member this.AminoAcids with get () =
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

    member this.ShortProfileSequences with get () =
        [|
            "CACGTG"
            "CACGTG"
            "CACGTG"
            "CACGTG"
        |]
    
    member this.LongProfileSequences with get () =
        [|
            "GTGGCTGCACCACGTGTATGC"
            "ACATCGCATCACGTGACCAGT"
            "CCTCGCACGTGGTGGTACAGT"
            "CTCGTTAGGACCATCACGTGA"
        |]

    member this.ThyminFocusedProfileSequences with get () =
        [|
            "GTGGCTGCACCACGTGTATGCCACGTG"
            "ACATCGCATCACGTGACCAGTTAGTTG"
            "CCTCGCACGTGGTGGTACAGTCGTACG"
            "GCATAAAGGACCATCACGTGAAGCTGC"
            "TTTTTTTTTTTTTTTTTTTTTTTTTTT"
        |]

    member this.VeryLongProfileSequences with get () =
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

    member this.FrequencyCompositeVector with get () =
        this.ShortProfileSequences
        |> Array.map (fun item -> BioArray.ofNucleotideString item)
        |> Array.map(fun bioArray -> createFCVOf bioArray)
       
    member this.ProbabilityCompositeVector with get () =
        this.FrequencyCompositeVector
        |> Array.map(fun bioArray -> createPCVOf bioArray)

    member this.PositionProbabilityMatrix with get () =
        this.ShortProfileSequences
        |> Array.map (fun item -> BioArray.ofNucleotideString item)
        |> Array.map(fun bioArray -> createPFMOf bioArray)
        |> Array.map(fun pfm -> createPPMOf pfm)
        |> Array.map(fun pbm -> normalizePPM this.ShortProfileSequences[0].Length this.DNABases this.Pseudocount pbm)

    member this.PositionWeightMatrix with get () =
        new PositionWeightMatrix(Array2D.length2 this.PositionProbabilityMatrix[0].Matrix)
