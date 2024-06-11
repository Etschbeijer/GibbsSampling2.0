namespace SiteSampler

open System
open FSharpAux
open BioFSharp
open HelperFunctions
open CompositeVector.Types
open CompositeVector.Functions
open PositionMatrix.Types
open PositionMatrix.Functions

module Functions =

    /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
    let getBestPWMSsWithBPV (motiveLength:int) (alphabet:#IBioItem[]) (source:BioArray.BioArray<#IBioItem>) (pcv:ProbabilityCompositeVector) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) (highValue:float) (highIndex:int) =
            if n + motiveLength > source.Length then log2(highValue), highIndex
            else
                let tmp =
                    let segment =
                        Array.skip n source
                        |> Array.take motiveLength
                    let pwMatrix = createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                    segment
                    |> calculateSegmentScoreBy pwMatrix
                if tmp > highValue then 
                    loop (n + 1) tmp n
                else 
                    loop (n + 1) highValue highIndex
        loop 0 0. 0

    /// Checks whether downstream of given positions a higher InformationContent is present or not. 
    /// If yes, the new InformationContent and positions are given back, otherwise the old ones.
    let getRightShiftedBestPWMSsWithBPV (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (startPositions:(float*int)[]) =
        let rnd = new Random()
        let source = [|0..sources.Length-1|]
        let randomSourceNumber = Array.shuffleFisherYates(rnd)(source)
        let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then 
                    acc
                else 
                    loop 0 acc (Array.copy acc)
            else                
                let unChosenArrays =
                    mergeArrays sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let unChosenStartPositions =  
                    mergeArrays bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                    |> Array.map2 (fun (source:BioArray.BioArray<#IBioItem>) (_, position) -> if position <= source.Length - motiveLength - 1 then position + 1 else position) unChosenArrays
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (getSegment motiveLength subSequence position) 
                         |> createPFMOf) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motiveLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else 
                        acc                    
                    )
                    bestmotive
        loop 0 startPositions startPositions

    /// Checks whether upstream of given positions a higher PositionWeightMatrixScore is present or not. 
    /// If yes, the new PositionWeightMatrixScore and positions are given back, otherwise the old ones.
    let getLeftShiftedBestPWMSsWithBPV (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (startPositions:(float*int)[]) =
        let rnd = new Random()
        let randomSourceNumber = Array.shuffleFisherYates(rnd)([|0..sources.Length-1|])
        let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
            if n = sources.Length then
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else 
                    loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =                   
                    mergeArrays bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                    |> Array.map (fun (_, position) -> if position > 0 then position - 1 else position)
                let unChosenArrays = 
                    mergeArrays sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (getSegment motiveLength subSequence position) 
                         |> createPFMOf) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motiveLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else 
                        acc                    
                    )
                    bestmotive
        loop 0 startPositions startPositions

    /// Checks the given Sequence for the existence of a conserved motive, by scoring each segment based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence.
    let findBestmotiveWithStartPosition (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (startPositions:(float*int)[]) =        
        let rnd = new Random()
        let randomSourceNumber = Array.shuffleFisherYates(rnd)([|0..sources.Length-1|])
        let rec loop (n:int) acc bestmotive =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else 
                    loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    mergeArrays acc.[0..randomSourceNumber.[n]-1] acc.[randomSourceNumber.[n]+1..]
                    |> Array.map (fun (_, position) -> position)
                let unChosenArrays =
                    mergeArrays sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (getSegment motiveLength subSequence position) 
                         |> createPFMOf) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motiveLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else 
                        acc                    
                    )
                    bestmotive
        loop 0 startPositions startPositions

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getPWMOfRandomStartsWithBPV (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =    
        let rnd = new Random()
        let randomSourceNumber = Array.shuffleFisherYates(rnd)([|0..sources.Length-1|])
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> List.toArray
            else
                let unChosenArrays =
                    mergeArrays (sources.[0..randomSourceNumber.[n]-1]) (sources.[randomSourceNumber.[n]+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        getRandomNumberInSequence motiveLength unChosen)
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (getSegment motiveLength subSequence position) 
                         |> createPFMOf) unChosenArrays randomStartPositions
                    |> fusePositionFrequencyMatrices motiveLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                loop (n + 1) (getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix::acc)
        loop 0 []

    /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
    /// Each iteration it is checked if better InformationContent is found or not.
    let getmotifsWithBestInformationContentWithBPV (numberOfRepetitions:int) (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =
        let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getPWMOfRandomStartsWithBPV motiveLength pseudoCount alphabet sources pcv
                            |> findBestmotiveWithStartPosition motiveLength pseudoCount alphabet sources pcv
                            |> getLeftShiftedBestPWMSsWithBPV motiveLength pseudoCount alphabet sources pcv
                            |> getRightShiftedBestPWMSsWithBPV motiveLength pseudoCount alphabet sources pcv
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|0., 0|]

    /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
    let getBestPWMSs (motiveLength:int) (alphabet:#IBioItem[]) (pseudoCount:float) (source:BioArray.BioArray<#IBioItem>) (fcVector:FrequencyCompositeVector) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) (highValue:float) (highIndex:int) =
            if n + motiveLength > source.Length then log2(highValue), highIndex
            else
                let tmp =
                    let segment =
                        Array.skip n source
                        |> Array.take motiveLength
                    let pcv =
                        increaseInPlaceFCVOf source fcVector
                        |> substractSegmentCountsFrom segment 
                        |> createNormalizedPCVOfFCV alphabet pseudoCount
                    let pwMatrix = createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                    segment
                    |> calculateSegmentScoreBy pwMatrix
                if tmp > highValue then loop (n + 1) tmp n
                else 
                    loop (n + 1) highValue highIndex
        loop 0 0. 0

    /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
    let getBestPWMSsWithPWM motiveLength (source:BioArray.BioArray<#IBioItem>) (pwMatrix:PositionWeightMatrix) =
        let rec loop (n:int) (highValue:float) (highIndex:int) =
            if n + motiveLength > source.Length then log2(highValue), highIndex
            else
                let tmp =
                    let segment =
                        Array.skip n source
                        |> Array.take motiveLength
                    segment
                    |> calculateSegmentScoreBy pwMatrix
                if tmp > highValue then loop (n + 1) tmp n
                else 
                    loop (n + 1) highValue highIndex
        loop 0 0. 0

    /// Checks whether downstream of given positions a higher InformationContent is present or not. 
    /// If yes, the new InformationContent and positions are given back, otherwise the old ones.
    let getRightShiftedBestPWMSs (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =
        let rnd = new Random()
        let randomSourceNumber = Array.shuffleFisherYates(rnd)([|0..sources.Length-1|])
        let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else 
                    loop 0 acc (Array.copy acc)
            else
                let unChosenArrays =
                    mergeArrays sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let unChosenStartPositions =
                    mergeArrays bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                    |> Array.map2 (fun (source:BioArray.BioArray<#IBioItem>) (_, position) -> if position <= source.Length - motiveLength - 1 then position + 1 else position) unChosenArrays
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) unChosenArrays unChosenStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        getSegment motiveLength subSequence position
                        |> createPFMOf) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motiveLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else 
                        acc                    
                    )
                    bestmotive
        loop 0 startPositions startPositions

    /// Checks whether upstream of given positions a higher PositionWeightMatrixScore is present or not. 
    /// If yes, the new PositionWeightMatrixScore and positions are given back, otherwise the old ones.
    let getLeftShiftedBestPWMSs (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =
        let rnd = new Random()
        let randomSourceNumber = Array.shuffleFisherYates(rnd)([|0..sources.Length-1|])
        let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
            if n = sources.Length then
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else 
                    loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    mergeArrays bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                    |> Array.map (fun (_, position) -> if position > 0 then position - 1 else position)
                let unChosenArrays =
                    mergeArrays sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) unChosenArrays unChosenStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        getSegment motiveLength subSequence position
                        |> createPFMOf) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motiveLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else 
                        acc                    
                    )
                    bestmotive
        loop 0 startPositions startPositions

    /// Checks the given sequence for the existence of a conserved motive, by scoring each segment based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated at each iteration if segments with higher scores are found until convergence.
    let getBestPWMSsWithStartPositions (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =        
        let rnd = new Random()
        let randomSourceNumber = Array.shuffleFisherYates(rnd)([|0..sources.Length-1|])
        let rec loop (n:int) acc bestmotive =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else 
                    loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    mergeArrays acc.[0..randomSourceNumber.[n]-1] acc.[randomSourceNumber.[n]+1..]
                    |> Array.map (fun (_, position) -> position)
                let unChosenArrays =
                    mergeArrays sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) unChosenArrays unChosenStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        getSegment motiveLength subSequence position
                        |> createPFMOf) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motiveLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else 
                        acc                    
                    )
                    bestmotive
        loop 0 startPositions startPositions

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getPWMOfRandomStarts (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        let rnd = new Random()
        let randomSourceNumber = Array.shuffleFisherYates(rnd)([|0..sources.Length-1|])
        let rec loop (n:int) (acc:((float*int)[])) =
            if n = sources.Length then acc
            else
                let unChosenArrays =
                    mergeArrays (sources.[0..randomSourceNumber.[n]-1]) (sources.[randomSourceNumber.[n]+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        getRandomNumberInSequence motiveLength unChosen)
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) unChosenArrays randomStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        getSegment motiveLength subSequence position
                        |> createPFMOf) unChosenArrays randomStartPositions
                    |> fusePositionFrequencyMatrices motiveLength
                    |> createPPMOf
                    |> normalizePPM (sources.Length - 1) alphabet pseudoCount
                loop (n + 1) (acc.[randomSourceNumber.[n]] <- getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                              acc)
        loop 0 (Array.zeroCreate sources.Length)

    /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
    /// Each iteration it is checked if better InformationContent is found or not.
    let getmotifsWithBestInformationContent (numberOfRepetitions:int) (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getPWMOfRandomStarts motiveLength pseudoCount alphabet sources
                            |> getBestPWMSsWithStartPositions motiveLength pseudoCount alphabet sources
                            |> getLeftShiftedBestPWMSs motiveLength pseudoCount alphabet sources
                            |> getRightShiftedBestPWMSs motiveLength pseudoCount alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|0., 0|]

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getmotifsWithBestPWMSOfPPM (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rnd = new Random()
        let randomSourceNumber = Array.shuffleFisherYates(rnd)([|0..sources.Length-1|])
        let rec loop (n:int) acc:((float*int)[]) =
            if n = sources.Length then acc
            else
                let unChosenArrays =
                    mergeArrays (sources.[0..randomSourceNumber.[n]-1]) (sources.[randomSourceNumber.[n]+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        getRandomNumberInSequence motiveLength unChosen)
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) 
                            unChosenArrays randomStartPositions
                    |> fuseFrequencyVectors alphabet
                loop (n + 1) (acc.[randomSourceNumber.[n]] <- getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                              acc)
        loop 0 (Array.zeroCreate sources.Length)

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getmotifsWithBestPWMSOfPWM (motiveLength:int) (sources:BioArray.BioArray<#IBioItem>[]) (pwm:PositionWeightMatrix) =
        let rnd = new Random()
        let randomSourceNumber = Array.shuffleFisherYates(rnd)([|0..sources.Length-1|])
        let rec loop (n:int) acc:((float*int)[]) =
            if n = sources.Length then acc
            else
                loop (n + 1) (acc.[randomSourceNumber.[n]] <- getBestPWMSsWithPWM motiveLength sources.[randomSourceNumber.[n]] pwm
                              acc)
        loop 0 (Array.zeroCreate sources.Length)

    /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
    /// Each iteration it is checked if better InformationContent is found or not.
    let getBestInformationContentOfPPM (numberOfRepetitions:int) (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getmotifsWithBestPWMSOfPPM motiveLength pseudoCount alphabet sources positionProbabilityMatrix
                            |> getBestPWMSsWithStartPositions motiveLength pseudoCount alphabet sources
                            |> getLeftShiftedBestPWMSs motiveLength pseudoCount alphabet sources
                            |> getRightShiftedBestPWMSs motiveLength pseudoCount alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|0., 0|]

    let doSiteSampling (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        getPWMOfRandomStarts motiveLength pseudoCount alphabet sources
        |> getBestPWMSsWithStartPositions motiveLength pseudoCount alphabet sources
        |> getLeftShiftedBestPWMSs motiveLength pseudoCount alphabet sources
        |> getRightShiftedBestPWMSs motiveLength pseudoCount alphabet sources

    let doSiteSamplingWithPCV (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =
        getPWMOfRandomStartsWithBPV motiveLength pseudoCount alphabet sources pcv
        |> findBestmotiveWithStartPosition motiveLength pseudoCount alphabet sources pcv
        |> getLeftShiftedBestPWMSsWithBPV motiveLength pseudoCount alphabet sources pcv
        |> getRightShiftedBestPWMSsWithBPV motiveLength pseudoCount alphabet sources pcv

    let doSiteSamplingWithPPM (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (ppM:PositionProbabilityMatrix) =
        getmotifsWithBestPWMSOfPPM motiveLength pseudoCount alphabet sources ppM
        |> getBestPWMSsWithStartPositions motiveLength pseudoCount alphabet sources
        |> getLeftShiftedBestPWMSs motiveLength pseudoCount alphabet sources
        |> getRightShiftedBestPWMSs motiveLength pseudoCount alphabet sources

    let doSiteSamplingWithPWM (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pwm:PositionWeightMatrix) =
        getmotifsWithBestPWMSOfPWM motiveLength sources pwm
        |> getBestPWMSsWithStartPositions motiveLength pseudoCount alphabet sources
        |> getLeftShiftedBestPWMSs motiveLength pseudoCount alphabet sources
        |> getRightShiftedBestPWMSs motiveLength pseudoCount alphabet sources

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getMotifsWithBestPWMSOfPPM (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> List.toArray
            else
                let unChosenArrays =
                    mergeArrays (sources.[0..n-1]) (sources.[n+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        getRandomNumberInSequence motifLength unChosen)
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motifLength position unchosenArray) 
                            unChosenArrays randomStartPositions
                    |> fuseFrequencyVectors alphabet
                loop (n + 1) (getBestPWMSs motifLength alphabet pseudoCount sources.[n] frequencyCompositeVector positionProbabilityMatrix::acc)
        loop 0 []