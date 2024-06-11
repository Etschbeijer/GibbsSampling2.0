namespace SuccessfullTests

open Xunit
open FluentAssertions
open BioFSharp
open HelperFunctions

open Fixture

module CompositeVectorTests =

    open CompositeVector.Functions
    
    [<Fact>]
    let ShouldCreateFrequenceCompositeVector () =
        
        let result =             
            BioArray.ofNucleotideString Fixture.ShortProfileSequences
            |> createFCVOf

        (result[Nucleotides.Nucleotide.A]).Should().Be(1, null, null)
        (result[Nucleotides.Nucleotide.T]).Should().Be(1, null, null)
        (result[Nucleotides.Nucleotide.G]).Should().Be(2, null, null)
        (result[Nucleotides.Nucleotide.C]).Should().Be(2, null, null)
        (result[Nucleotides.Nucleotide.Gap]).Should().Be(0, null, null)

    [<Fact>]
    let ShouldCreateProbabilityCompositeVector () =
        
        let result =             
            BioArray.ofNucleotideString Fixture.ShortProfileSequences
            |> createFCVOf
            |> createPCVOf

        (result[Nucleotides.Nucleotide.A]).Should().Be(1., null, null)
        (result[Nucleotides.Nucleotide.T]).Should().Be(1., null, null)
        (result[Nucleotides.Nucleotide.G]).Should().Be(2., null, null)
        (result[Nucleotides.Nucleotide.C]).Should().Be(2., null, null)
        (result[Nucleotides.Nucleotide.Gap]).Should().Be(0., null, null)

    [<Fact>]
    let ShouldCreateNormalizedProbabilityCompositeVector () =
        
        let fcv =             
            BioArray.ofNucleotideString Fixture.ShortProfileSequences
            |> createFCVOf

        let pcv =  createPCVOf fcv

        let result = createNormalizedPCVOfFCV Fixture.DNABases Fixture.Pseudocount fcv

        let calcNormalizaion(nultetide: Nucleotides.Nucleotide) = 
            (pcv[nultetide] + Fixture.Pseudocount) / float ((Array.sum fcv.Array) + Fixture.DNABases.Length)

        (result[Nucleotides.Nucleotide.A]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.A, null, null)
        (result[Nucleotides.Nucleotide.T]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.T, null, null)
        (result[Nucleotides.Nucleotide.G]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.G, null, null)
        (result[Nucleotides.Nucleotide.C]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.C, null, null)
        (result[Nucleotides.Nucleotide.Gap]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.Gap, null, null)

    [<Fact>]
    let ShouldCalculateSegmentScore () =
        
        let ba =            
            BioArray.ofNucleotideString Fixture.ShortProfileSequences
            |> (fun item -> item[0..2])

        let fcv =             
            BioArray.ofNucleotideString Fixture.ShortProfileSequences
            |> createFCVOf

        let pncv =             
            createNormalizedPCVOfFCV Fixture.DNABases Fixture.Pseudocount fcv

        let result = calculateSegmentScoreBy pncv ba

        let calcSegmentValue =
            let pcv = pncv
            pcv[ba[0]] * pcv[ba[1]] * pcv[ba[2]]

        result.Should().Be(calcSegmentValue, null, null)



module PositionMatrixTests =

    open CompositeVector.Functions
    open PositionMatrix.Functions
    
    [<Fact>]
    let ShouldCreatePositionFrequencyMatrix () =
    
        let fixture = new Fixture()
        
        let result =             
            BioArray.ofNucleotideString Fixture.ShortProfileSequences
            |> createPFMOf

        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 0].Should().Be(1, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.A |> int) - 42, 1].Should().Be(1, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 2].Should().Be(1, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 3].Should().Be(1, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.T |> int) - 42, 4].Should().Be(1, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 5].Should().Be(1, null, null)

    [<Fact>]
    let ShouldCreatePositionProbabilityMatrix () =
    
        let fixture = new Fixture()
        
        let result =             
            BioArray.ofNucleotideString Fixture.ShortProfileSequences
            |> createPFMOf
            |> createPPMOf

        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 0].Should().Be(1., null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.A |> int) - 42, 1].Should().Be(1., null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 2].Should().Be(1., null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 3].Should().Be(1., null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.T |> int) - 42, 4].Should().Be(1., null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 5].Should().Be(1., null, null)

    [<Fact>]
    let ShouldNormalizePositionProbabilityMatrix () =
        
        let fcv =             
            BioArray.ofNucleotideString Fixture.ShortProfileSequences
            |> createFCVOf

        let pcv = createPCVOf fcv

        let result =             
            BioArray.ofNucleotideString Fixture.ShortProfileSequences
            |> createPFMOf
            |> createPPMOf
            |> normalizePPM Fixture.ShortProfileSequences.Length Fixture.DNABases Fixture.Pseudocount

        let calcNormalizaion(nultetide: Nucleotides.Nucleotide) = 
            let value = 
                if pcv[nultetide] > 1 then 
                    1.
                else 
                    pcv[nultetide]
            (value + Fixture.Pseudocount) / float ((Array.sum fcv.Array) + Fixture.DNABases.Length)

        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 0].Should().Be(calcNormalizaion Nucleotides.Nucleotide.C, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.A |> int) - 42, 1].Should().Be(calcNormalizaion Nucleotides.Nucleotide.A, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 2].Should().Be(calcNormalizaion Nucleotides.Nucleotide.C, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 3].Should().Be(calcNormalizaion Nucleotides.Nucleotide.G, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.T |> int) - 42, 4].Should().Be(calcNormalizaion Nucleotides.Nucleotide.T, null, null)
        result.Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 5].Should().Be(calcNormalizaion Nucleotides.Nucleotide.G, null, null)

module SiteSamplerTests =

    open SiteSampler.Functions

    [<Fact>]
    let ShouldGetBestSegmet () =
        
        let source = BioArray.ofNucleotideString Fixture.ShortProfileSequences

        let result = 
            getBestPWMSsWithBPV 3 Fixture.DNABases source Fixture.ProbabilityCompositeVector Fixture.PositionProbabilityMatrix

        (snd result).Should().Be(0, null, null)

    [<Fact>]
    let ShouldGetNotGetRightSegmet () =

        let fixture = new Fixture()
        
        let source = BioArray.ofNucleotideString Fixture.ShortProfileSequences

        let result = 
            getRightShiftedBestPWMSsWithBPV 3 1 Fixture.DNABases [|source; source; source; source|] Fixture.ProbabilityCompositeVector [|(5., 0); (5., 1); (5., 2); (5., 0)|]

        (snd result[0]).Should().Be(0, null, null)

    [<Fact>]
    let ShouldGetNotGetLeftSegmet () =
        
        let source = BioArray.ofNucleotideString Fixture.ShortProfileSequences

        let result = 
            getLeftShiftedBestPWMSsWithBPV 3 1 Fixture.DNABases [|source; source; source; source|] Fixture.ProbabilityCompositeVector [|(5., 0); (5., 1); (5., 2); (5., 0)|]

        (snd result[0]).Should().Be(0, null, null)

    [<Fact>]
    let ShouldGetMotifsWithBestPWM () =
        
        let source = BioArray.ofNucleotideString Fixture.ShortProfileSequences

        let result = 
            getMotifsWithBestPWMSOfPPM 3 1 Fixture.DNABases [|source; source; source; source|] Fixture.PositionProbabilityMatrix

        (snd result[0]).Should().Be(0, null, null)

module MotifSamplerTests =

    open MotifSampler.Functions

    [<Fact>]
    let ShouldCalculatePWMsForSegmentCombinations () =
        
        let result = calculatePWMsForSegmentCombinations 1 3 1 [(2, 0); (3, 1); (4, 0); (3, 2)]

        (System.Math.Round(result[0].PWMS, 2)).Should().Be(1.58)
        result[0].Positions[0].Should().Be(1)
        result[1].PWMS.Should().Be(2)
        result[1].Positions[0].Should().Be(0)
        (System.Math.Round(result[2].PWMS, 2)).Should().Be(1.58)
        result[2].Positions[0].Should().Be(2)

    [<Fact>]
    let ShouldCalculateNormalizedSegmentScores () =
        
        let result = calculateNormalizedSegmentScores 1 3 3 Fixture.DNABases Fixture.ProbabilityCompositeVector Fixture.PositionWeightMatrix

        result[0].PWMS.Should().Be(2)
        result[1].PWMS.Should().Be(4)
        result[2].PWMS.Should().Be(0)

module HeuristicTests =

    open SiteSampler.Functions

    /// Heuristic function, smae result is possible but should not happen often
    [<Fact>]
    let ShouldGetDifferentSegments () =
        
        let fixture = new Fixture()

        let complexProfileSequences =
            Fixture.ComplexProfileSequences
            |> Array.map (fun bioArray -> BioArray.ofNucleotideString bioArray)

        let resultI =             
            Array.init 10 (fun _ -> getmotifsWithBestInformationContent 10 6 0.0001 Fixture.DNABases complexProfileSequences)
        
        let orderedResultI =
            resultI
            |> Array.countBy (fun items -> items |> Array.map (fun item -> snd item))
            |> Array.sortByDescending (fun (_, i) -> i)

        let resultII =             
            Array.init 10 (fun _ -> getmotifsWithBestInformationContent 10 6 0.0001 Fixture.DNABases complexProfileSequences)
        
        let orderedResultII =
            resultII
            |> Array.countBy (fun items -> items |> Array.map (fun item -> snd item))
            |> Array.sortByDescending (fun (_, i) -> i)

        orderedResultI.Should().NotBeSameAs(orderedResultII)

    [<Fact>]
    let ShouldMergeArrays () =

        let arrayI = [|1; 2; 3;|]
        let arrayII = [|4; 4; 4; 5; 6|]         

        let mergedArrays = mergeArrays arrayI arrayII

        mergedArrays
