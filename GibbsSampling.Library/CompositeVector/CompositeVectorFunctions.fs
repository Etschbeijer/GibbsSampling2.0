namespace CompositeVector

open BioFSharp
open CompositeVector.Types

module Functions =

    /// Increase counter of position by 1.
    let increaseInPlaceFCV (bioItem:#IBioItem) (frequencyCompositeVector:FrequencyCompositeVector) = 
        frequencyCompositeVector.[bioItem] <- frequencyCompositeVector.[bioItem] + 1
        frequencyCompositeVector

    /// Increase counter of position by n.
    let increaseInPlaceFCVBy (bioItem:'a) (frequencyCompositeVector:FrequencyCompositeVector) n = 
        frequencyCompositeVector.[bioItem] <- frequencyCompositeVector.[bioItem] + n
        frequencyCompositeVector

    /// Create a FrequencyCompositeVector based on BioArrays and exclude the specified segments.
    let createFCVOf (resSources:BioArray.BioArray<#IBioItem>) =
        let backGroundCounts = new FrequencyCompositeVector()   
        Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts resSources

    /// Create new FrequencyCompositeVector based on the sum of the positions of an array of FrequencyVectors.
    let fuseFrequencyVectors (alphabet:#IBioItem[]) (bfVectors:FrequencyCompositeVector[]) =
        let backgroundFrequencyVector = new FrequencyCompositeVector()
        for fcVector in bfVectors do
            for bioItem in alphabet do
                backgroundFrequencyVector.[bioItem] <- backgroundFrequencyVector.[bioItem] + (fcVector.[bioItem])
        backgroundFrequencyVector

    /// Create FrequencyCompositeVector based on BioArray and excludes the specified segment.
    let createFCVWithout (motiveLength:int) (position:int) (resSource:BioArray.BioArray<#IBioItem>) =
        let backGroundCounts = new FrequencyCompositeVector()
        Array.append resSource.[0..(position-1)] resSource.[(position+motiveLength)..]
        |> Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts

    /// Create FrequencyCompositeVector based on BioArray.
    let increaseInPlaceFCVOf (resSources:BioArray.BioArray<#IBioItem>) (backGroundCounts:FrequencyCompositeVector) =   
        resSources
        |> Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts

    /// Subtracts the amount of elements in the given source from the FrequencyCompositeVector.
    let substractSegmentCountsFrom (source:BioArray.BioArray<#IBioItem>) (fcVector:FrequencyCompositeVector) =
        let bfVec = new FrequencyCompositeVector(fcVector.Array)
        for bioItem in source do
            bfVec.[bioItem] <- (if fcVector.[bioItem] - 1 > 0 then fcVector.[bioItem] - 1 else 0)
        bfVec

    /// Increase counter of position by 1.
    let increaseInPlacePCV (bioItem:'a) (backGroundProbabilityVector:ProbabilityCompositeVector) = 
        backGroundProbabilityVector.[bioItem] <- backGroundProbabilityVector.[bioItem] + 1.
        backGroundProbabilityVector

    /// Increase counter of position by n.
    let increaseInPlacePCVBy (bioItem:'a) n (backGroundProbabilityVector:ProbabilityCompositeVector) = 
        backGroundProbabilityVector.[bioItem] <- backGroundProbabilityVector.[bioItem] + n
        backGroundProbabilityVector

    /// Create a ProbabilityCompositeVector based on a FrequencyCompositeVector by replacing the integers with floats.
    let createPCVOf (caArray:FrequencyCompositeVector) =
        caArray.Array
        |> Array.map (fun item -> float item)
        |> fun item -> new ProbabilityCompositeVector(item)

    /// Create normalized ProbabilityCompositeVector based on FrequencyCompositeVector.
    let createNormalizedPCVOfFCV (alphabet:#IBioItem[]) (pseudoCount:float) (frequencyCompositeVector:FrequencyCompositeVector) =
        let backGroundProbabilityVector = createPCVOf frequencyCompositeVector
        let sum = (float (Array.sum frequencyCompositeVector.Array)) + ((float alphabet.Length) * pseudoCount)
        for item in alphabet do
            backGroundProbabilityVector.[item] <- (backGroundProbabilityVector.[item] + pseudoCount)/sum
        backGroundProbabilityVector

    /// Calculate the score of the given sequence based on the probabilities of the ProbabilityCompositeVector.
    let calculateSegmentScoreBy (pcv:ProbabilityCompositeVector) (bioItems:BioArray.BioArray<#IBioItem>) =
        Array.fold (fun (value:float) (bios:#IBioItem) -> value * (pcv.[bios])) 1. bioItems
