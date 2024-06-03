namespace MotifSampler

open FSharpAux
open BioFSharp
open HelperFunctions
open CompositeVector.Types
open CompositeVector.Functions
open PositionMatrix.Types
open PositionMatrix.Functions
open SiteSampler.Functions
open MotifSampler.Types

module Functions =
    
    /// Creates a MotifIndex based on the probability of the given segement(s) start position(s).
    let createMotifIndex pwms pos =
        {
            MotifIndex.PWMS         = pwms
            MotifIndex.Positions    = pos
        }

    /// Calculates all possible segment combination of m segments, which do not overlap in the given width 
    /// and gives you the probability for it.
    let calculatePWMsForSegmentCombinations (cutOff:float) (width:int) (m:int) (set:list<float*int>) =
        let rec loop prob positions size set = 
            seq
                {
                    match size, set with
                    | n, x::xs
                        when n >= 0 ->
                            if checkForDistance width (snd x::positions) then
                                if log2(fst x*prob) > cutOff then
                                    yield! loop (fst x*prob) (snd x::positions) (n - 1) xs
                            if n >= 0 then 
                                yield! loop prob positions n xs
                    | 0, [] -> yield createMotifIndex (log2(prob)) positions
                    | _, [] -> ()
                    | _ -> () 
                }
        loop 1. [] m set
        |> List.ofSeq

    /// Normalizes the probabilities of all MotifMemories to the sum of all probabilities and picks one by random, 
    /// but those with a higher value have a higher chance to get picked.
    let rouletteWheelSelection (pick:float) (items:MotifIndex list) =
        let normalizedItems =
            let sum = List.sum (items |> List.map (fun item -> item.PWMS))
            items
            |> List.map (fun item -> item.PWMS/sum, item)
        let rec loop acc n =
            if acc <= pick && pick <= acc + fst normalizedItems.[n] then 
                items.[n]
            else 
                loop (acc + fst normalizedItems.[n]) (n + 1)
        loop 0. 0       

    /// Calculates the normalized segment score based on the given PositionWeightMatrix and 
    /// BackGroundProbabilityVecor of all segment combinations of a given sequence. 
    /// The amount of combinations is given by the amount of motifs.
    let calculateNormalizedSegmentScores (cutOff:float) (motifAmount:int) (motifLength:int) (source:BioArray.BioArray<#IBioItem>) (pcv:ProbabilityCompositeVector) (pwMatrix:PositionWeightMatrix) =
        let segments =
            let rec loop n acc =
                if n + motifLength = source.Length+1 then 
                    List.rev acc
                else
                    let tmp =
                        source
                        |> Array.skip n
                        |> (fun items -> Array.take motifLength items, n)
                    loop (n+1) (tmp::acc)
            loop 0 []
        let segmentScores =
            segments
            |> List.map (fun segment -> 
                calculateSegmentScoreBy pwMatrix (fst segment), (snd segment))
        let backGroundScores =
            segments
            |> List.map (fun segment -> CompositeVector.Functions.calculateSegmentScoreBy pcv (fst segment))
            |> List.map (fun segmentScore -> createMotifIndex segmentScore [])
        let tmp =
            let rec loop n acc =
                if n > motifAmount then 
                    List.rev acc
                else loop (n+1) (calculatePWMsForSegmentCombinations cutOff motifLength n segmentScores::acc)
            loop 1 [backGroundScores]
        tmp
        |> List.concat

    /// Checks the given sequence for the existence of conserved motif(s), by scoring each segment and segment combination based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence. 
    let findBestMotifPositionsWithStartPositionByPCV (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (motifMem:MotifIndex[]) =
        let rec loop (n:int) acc bestMotif =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> item.Positions)) = (bestMotif |> Array.map (fun item -> item.Positions)) then 
                    acc
                else 
                    loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    mergeArrays acc.[0..n-1] acc.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    mergeArrays sources.[0..n-1] sources.[n+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            getSegment motifLength subSequence position
                            |> createPFMOf)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fusePositionFrequencyMatrices motifLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motifAmount motifLength sources.[n] pcv positionWeightMatrix
                    |> List.sortByDescending (fun item -> item.PWMS)
                    |> List.head
                loop 
                    (n + 1) 
                    (if tmp.PWMS > acc.[n].PWMS then 
                        acc.[n] <- tmp
                        acc
                     else 
                        acc                    
                    )
                    bestMotif
        loop 0 (Array.copy motifMem) (Array.copy motifMem)


    /// Checks the given Sequence for the existence of conserved motif(s), by scoring each segment and segment combination based on the given start positions.
    /// The nex segments for the nex PositionWeightMatrix are picked by chance after calculating the segment scores for each possible combination 
    /// but those with higher scores have a higher chance to be picked.
    let findBestMotifPositionsWithStartPositionsByPCV (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (motifMem:MotifIndex[]) =
        let rnd = new System.Random()
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> Array.ofList
            else
                let unChosenStartPositions =
                    mergeArrays motifMem.[0..n-1] motifMem.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    mergeArrays sources.[0..n-1] sources.[n+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            getSegment motifLength subSequence position
                            |> createPFMOf)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fusePositionFrequencyMatrices motifLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motifAmount motifLength sources.[n] pcv positionWeightMatrix
                    |> rouletteWheelSelection (rnd.NextDouble())
                loop (n+1) (tmp::acc)
        loop 0 []

    ///Repeats the motif sampler until convergence or until a given number of repetitions is done.
    let findBestInormationContentContainingMotifsWithPCV (numberOfRepetitions:int) (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =
        let rec loop (n:int) (acc:MotifIndex[]) (bestPWMS:MotifIndex[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getPWMOfRandomStartsWithBPV motifLength pseudoCount alphabet sources pcv
                            |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
                            |> findBestMotifPositionsWithStartPositionsByPCV motifAmount motifLength pseudoCount cutOff alphabet sources pcv
                            |> findBestMotifPositionsWithStartPositionByPCV motifAmount motifLength pseudoCount cutOff alphabet sources pcv
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|createMotifIndex 0. []|]

    /// Checks the given Sequence for the existence of conserved motif(s), by scoring each segment and segment combination based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence. 
    let findBestMotifIndicesWithStartPositions (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (motifMem:MotifIndex[]) =
        let rec loop (n:int) acc bestMotif =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> item.Positions)) = (bestMotif |> Array.map (fun item -> item.Positions)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    mergeArrays acc.[0..n-1] acc.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    mergeArrays sources.[0..n-1] sources.[n+1..]
                let backgroundProbabilityVector =
                    Array.map2 (fun array positions -> 
                        positions
                        |> Array.map (fun position -> 
                            createFCVWithout motifLength position array)
                                ) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fuseFrequencyVectors alphabet
                    |> increaseInPlaceFCVOf sources.[n]
                    |> createNormalizedPCVOfFCV alphabet pseudoCount
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            getSegment motifLength subSequence position
                            |> createPFMOf)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fusePositionFrequencyMatrices motifLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motifAmount motifLength sources.[n] backgroundProbabilityVector positionWeightMatrix
                    |> List.sortByDescending (fun item -> item.PWMS)
                    |> List.head
                loop 
                    (n + 1) 
                    (if tmp.PWMS > acc.[n].PWMS then 
                        acc.[n] <- tmp
                        acc
                     else 
                        acc                    
                    )
                    bestMotif
        loop 0 (Array.copy motifMem) (Array.copy motifMem)


    /// Checks the given Sequence for the existence of conserved motif(s), by scoring each segment and segment combination based on the given start positions.
    /// The nex segments for the nex PositionWeightMatrix are picked by chance after calculating the segment scores for each possible combination 
    /// but those with higher scores have a higher chance to be picked.
    let findBestMotifIndicesByWithStartPositions (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (motifMem:MotifIndex[]) =
        let rnd = new System.Random()
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> Array.ofList
            else
                let unChosenStartPositions =
                    mergeArrays motifMem.[0..n-1] motifMem.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    mergeArrays sources.[0..n-1] sources.[n+1..]
                let backgroundProbabilityVector =
                    Array.map2 (fun array positions -> 
                        positions
                        |> Array.map (fun position -> 
                            createFCVWithout motifLength position array)
                                ) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fuseFrequencyVectors alphabet
                    |> increaseInPlaceFCVOf sources.[n]
                    |> createNormalizedPCVOfFCV alphabet pseudoCount
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            getSegment motifLength subSequence position
                            |> createPFMOf)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fusePositionFrequencyMatrices motifLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motifAmount motifLength sources.[n] backgroundProbabilityVector positionWeightMatrix
                    |> rouletteWheelSelection (rnd.NextDouble())
                loop (n+1) (tmp::acc)
        loop 0 []

    ///Repeats the motif sampler until convergence or until a given number of repetitions is done.
    let getMotifsWithBestInformationContents (numberOfRepetitions:int) (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        let rec loop (n:int) (acc:MotifIndex[]) (bestPWMS:MotifIndex[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getPWMOfRandomStarts motifLength pseudoCount alphabet sources
                            |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
                            |> findBestMotifIndicesByWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources
                            |> findBestMotifIndicesWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|createMotifIndex 0. []|]

    ///Repeats the motif sampler until convergence or until a given number of repetitions is done.
    let getBestPWMSsOfPPM (numberOfRepetitions:int) (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) (acc:MotifIndex[]) (bestPWMS:MotifIndex[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getMotifsWithBestPWMSOfPPM motifLength pseudoCount alphabet sources positionProbabilityMatrix
                            |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
                            |> findBestMotifIndicesByWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources
                            |> findBestMotifIndicesWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|createMotifIndex 0. []|]

    let doMotifSamplingWithPPM (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        getMotifsWithBestPWMSOfPPM motifLength pseudoCount alphabet sources positionProbabilityMatrix
        |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
        |> findBestMotifIndicesByWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources
        |> findBestMotifIndicesWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources

    let doMotifSampling (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        getPWMOfRandomStarts motifLength pseudoCount alphabet sources
        |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
        |> findBestMotifIndicesByWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources
        |> findBestMotifIndicesWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources
