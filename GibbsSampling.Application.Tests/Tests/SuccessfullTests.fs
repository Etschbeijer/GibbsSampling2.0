namespace SuccessfullTests

open Xunit
open FluentAssertions
open BioFSharp

open Fixture

module CompositeVectorTests =

    open CompositeVector.Functions
    
    [<Fact>]
    let ShouldCreateFrequenceCompositeVector () =

        let fixture = new Fixture()
        
        let result = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)
            |> Array.map(fun bioArray -> createFCVOf bioArray)

        (result[0][Nucleotides.Nucleotide.A]).Should().Be(1, null, null)
        (result[0][Nucleotides.Nucleotide.T]).Should().Be(1, null, null)
        (result[0][Nucleotides.Nucleotide.G]).Should().Be(2, null, null)
        (result[0][Nucleotides.Nucleotide.C]).Should().Be(2, null, null)
        (result[0][Nucleotides.Nucleotide.Gap]).Should().Be(0, null, null)

    [<Fact>]
    let ShouldCreateProbabilityCompositeVector () =
    
        let fixture = new Fixture()
        
        let result = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)
            |> Array.map(fun bioArray -> createFCVOf bioArray)
            |> Array.map(fun bioArray -> createPCVOf bioArray)

        (result[0][Nucleotides.Nucleotide.A]).Should().Be(1., null, null)
        (result[0][Nucleotides.Nucleotide.T]).Should().Be(1., null, null)
        (result[0][Nucleotides.Nucleotide.G]).Should().Be(2., null, null)
        (result[0][Nucleotides.Nucleotide.C]).Should().Be(2., null, null)
        (result[0][Nucleotides.Nucleotide.Gap]).Should().Be(0., null, null)

    [<Fact>]
    let ShouldCreateNormalizedProbabilityCompositeVector () =
    
        let fixture = new Fixture()
        
        let fcv = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)
            |> Array.map(fun bioArray -> createFCVOf bioArray)

        let pcv = 
            fcv
            |> Array.map(fun bioArray -> createPCVOf bioArray)

        let result = 
            fcv
            |> Array.map(fun bioArray -> createNormalizedPCVOfFCV fixture.DNABases fixture.Pseudocount bioArray)

        let calcNormalizaion(nultetide: Nucleotides.Nucleotide) = 
            (pcv[0][nultetide] + fixture.Pseudocount) / float ((Array.sum fcv[0].Array) + fixture.DNABases.Length)

        (result[0][Nucleotides.Nucleotide.A]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.A, null, null)
        (result[0][Nucleotides.Nucleotide.T]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.T, null, null)
        (result[0][Nucleotides.Nucleotide.G]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.G, null, null)
        (result[0][Nucleotides.Nucleotide.C]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.C, null, null)
        (result[0][Nucleotides.Nucleotide.Gap]).Should().Be(calcNormalizaion Nucleotides.Nucleotide.Gap, null, null)

    [<Fact>]
    let ShouldCalculateSegmentScore () =
    
        let fixture = new Fixture()
        
        let ba =
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)
            |> Array.map (fun item -> item[0..2])

        let fcv = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)
            |> Array.map(fun bioArray -> createFCVOf bioArray)

        let pncv = 
            fcv
            |> Array.map(fun bioArray -> createNormalizedPCVOfFCV fixture.DNABases fixture.Pseudocount bioArray)

        let result =
            Array.map2(fun probabilityArray b -> calculateSegmentScoreBy probabilityArray b) pncv ba

        let calcSegmentValue =
            let pcv = pncv[0]
            pcv[ba[0][0]] * pcv[ba[0][1]] * pcv[ba[0][2]]

        result[0].Should().Be(calcSegmentValue, null, null)



module PositionMatrixTests =

    open CompositeVector.Functions
    open PositionMatrix.Functions
    
    [<Fact>]
    let ShouldCreatePositionFrequencyMatrix () =
    
        let fixture = new Fixture()
        
        let result = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)
            |> Array.map(fun bioArray -> createPFMOf bioArray)

        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 0].Should().Be(1, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.A |> int) - 42, 1].Should().Be(1, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 2].Should().Be(1, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 3].Should().Be(1, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.T |> int) - 42, 4].Should().Be(1, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 5].Should().Be(1, null, null)

    [<Fact>]
    let ShouldCreatePositionProbabilityMatrix () =
    
        let fixture = new Fixture()
        
        let result = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)
            |> Array.map(fun bioArray -> createPFMOf bioArray)
            |> Array.map(fun pfm -> createPPMOf pfm)

        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 0].Should().Be(1., null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.A |> int) - 42, 1].Should().Be(1., null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 2].Should().Be(1., null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 3].Should().Be(1., null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.T |> int) - 42, 4].Should().Be(1., null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 5].Should().Be(1., null, null)

    [<Fact>]
    let ShouldNormalizePositionProbabilityMatrix () =
    
        let fixture = new Fixture()
        
        let fcv = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)
            |> Array.map(fun bioArray -> createFCVOf bioArray)

        let pcv = 
            fcv
            |> Array.map(fun bioArray -> createPCVOf bioArray)

        let result = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)
            |> Array.map(fun bioArray -> createPFMOf bioArray)
            |> Array.map(fun pfm -> createPPMOf pfm)
            |> Array.map(fun pbm -> normalizePPM fixture.ShortProfileSequences[0].Length fixture.DNABases fixture.Pseudocount pbm)

        let calcNormalizaion(nultetide: Nucleotides.Nucleotide) = 
            let value = 
                if pcv[0][nultetide] > 1 then 1.
                else pcv[0][nultetide]
            (value + fixture.Pseudocount) / float ((Array.sum fcv[0].Array) + fixture.DNABases.Length)

        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 0].Should().Be(calcNormalizaion Nucleotides.Nucleotide.C, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.A |> int) - 42, 1].Should().Be(calcNormalizaion Nucleotides.Nucleotide.A, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.C |> int) - 42, 2].Should().Be(calcNormalizaion Nucleotides.Nucleotide.C, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 3].Should().Be(calcNormalizaion Nucleotides.Nucleotide.G, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.T |> int) - 42, 4].Should().Be(calcNormalizaion Nucleotides.Nucleotide.T, null, null)
        result[0].Matrix[(BioItem.symbol Nucleotides.Nucleotide.G |> int) - 42, 5].Should().Be(calcNormalizaion Nucleotides.Nucleotide.G, null, null)

module SiteSamplerTests =

    open SiteSampler.Functions

    [<Fact>]
    let ShouldGetBestSegmet () =

        let fixture = new Fixture()
        
        let source = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)

        let result = 
            getBestPWMSsWithBPV 3 fixture.DNABases source[0] fixture.ProbabilityCompositeVector[0] fixture.PositionProbabilityMatrix[0]

        (snd result).Should().Be(0, null, null)

    [<Fact>]
    let ShouldGetNotGetRightSegmet () =

        let fixture = new Fixture()
        
        let source = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)

        let result = 
            getRightShiftedBestPWMSsWithBPV 3 1 fixture.DNABases source fixture.ProbabilityCompositeVector[0] [|(5., 0); (5., 1); (5., 2); (5., 0)|]

        (snd result[0]).Should().Be(0, null, null)

    [<Fact>]
    let ShouldGetNotGetLeftSegmet () =

        let fixture = new Fixture()
        
        let source = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)

        let result = 
            getLeftShiftedBestPWMSsWithBPV 3 1 fixture.DNABases source fixture.ProbabilityCompositeVector[0] [|(5., 0); (5., 1); (5., 2); (5., 0)|]

        (snd result[0]).Should().Be(0, null, null)

    [<Fact>]
    let ShouldGetMotifsWithBestPWM () =

        let fixture = new Fixture()
        
        let source = 
            fixture.ShortProfileSequences
            |> Array.map (fun item -> BioArray.ofNucleotideString item)

        let result = 
            getMotifsWithBestPWMSOfPPM 3 1 fixture.DNABases source fixture.PositionProbabilityMatrix[0]

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

        let fixture = new Fixture()
        
        let result = calculateNormalizedSegmentScores 1 3 3 fixture.DNABases fixture.ProbabilityCompositeVector[0] fixture.PositionWeightMatrix

        result[0].PWMS.Should().Be(2)
        result[1].PWMS.Should().Be(4)
        result[2].PWMS.Should().Be(0)
