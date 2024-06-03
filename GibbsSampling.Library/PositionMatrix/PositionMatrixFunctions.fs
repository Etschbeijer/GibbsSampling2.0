namespace PositionMatrix

open FSharpAux
open BioFSharp
open CompositeVector.Types
open PositionMatrix.Types

module Functions =

    ///Checks whether all elements in the list have a wider distance than width or not.
    let internal checkForDistance (width:int) (items:int list) =
        if items.Length <= 0 || items.Length = 1 then true
        else
            let rec loop n i =
                if n = items.Length-1 then true
                else
                    if i >= items.Length then loop (n+1) (n+2)
                    else
                        if Operators.abs(items.[n]-items.[i]) > width then
                            loop n (i+1)
                        else 
                            false
            loop 0 1

    /// Get an integer which is between 0 and the length of the sequence - segmentLength
    let internal getRandomNumberInSequence (segmentLength:int) (source:'a[]) =
        let rnd = System.Random()
        rnd.Next(0, source.Length-segmentLength+1)

    /// Create a specific sub sequence of the source sequence based on the given length and starting position. Do not forget, sequences start counting at 0!
    let internal getSegment (subsequenceLength:int) (source:'a[]) (startPoint:int) =
        source 
        |> Array.skip startPoint 
        |> Array.take subsequenceLength 
        |> Array.ofSeq

    ///Finds the best information content in an array of arrays of PWMSs and positions.
    let getBestInformationContent (item:((float*int)[])[]) =
        let rec loop (n:int) (bestPWMS:(float*int)[]) =
            if n = item.Length then bestPWMS
            else
                let informationContentItem =
                    item.[n]
                    |> Array.fold (fun baseValue (pwms, _) -> pwms + baseValue) 0.
                let informationContentBestPWMS =
                    bestPWMS
                    |> Array.fold (fun baseValue (pwms, _) -> pwms + baseValue) 0.
                if informationContentItem > informationContentBestPWMS then
                    loop (n + 1) item.[n]
                else
                    loop (n + 1) bestPWMS
        loop 0 [||]

    /// Increase counter of PositionFrequencyMatrix at fixed position by 1.
    let increaseInPlacePFM (pos:int) (bioItem:'a when 'a :> IBioItem) (positionFrequencyMatrix:PositionFrequencyMatrix) = 
        positionFrequencyMatrix.[bioItem, pos] <- positionFrequencyMatrix.[bioItem, pos] + 1
        positionFrequencyMatrix

    /// Increase counter of PositionFrequencyMatrix at fixed position by n.
    let increaseInPlacePFMBy (pos:int) (bioItem:'a when 'a :> IBioItem) n (positionFrequencyMatrix:PositionFrequencyMatrix) = 
        positionFrequencyMatrix.[bioItem, pos] <- positionFrequencyMatrix.[bioItem, pos] + n
        positionFrequencyMatrix

    /// Create PositionFrequencyMatrix based on BioArray.
    let createPFMOf (source:BioArray.BioArray<#IBioItem>) =
        let positionFrequencyMatrix = new PositionFrequencyMatrix(source.Length)
        source
        |> Array.foldi (fun row (cm) column -> increaseInPlacePFM row column cm) (positionFrequencyMatrix)

    /// Create new PositionFrequencyMatrix based on the sum of the positions of an array of Countmatrices.
    let fusePositionFrequencyMatrices (motiveLength:int) (countMatrices:PositionFrequencyMatrix[]) =
        let positionFrequencyMatrix = new PositionFrequencyMatrix(motiveLength)
        if countMatrices.Length > 0 then
            for cMatrix in countMatrices do
                for column = 0 to (Array2D.length1 cMatrix.Matrix)-1 do
                        for row = 0 to (Array2D.length2 cMatrix.Matrix)-1 do
                            positionFrequencyMatrix.Matrix.[column, row] <- positionFrequencyMatrix.Matrix.[column, row] + (cMatrix.Matrix.[column, row])
            positionFrequencyMatrix
        else 
            positionFrequencyMatrix

    /// Increase counter of PositionProbabilityMatrix at fixed position by 1.
    let increaseInPlacePPM (pos:int) (bioItem:'a when 'a :> IBioItem) (positionProbabilityMatrix:PositionProbabilityMatrix) = 
        positionProbabilityMatrix.[bioItem, pos] <- positionProbabilityMatrix.[bioItem, pos] + 1.
        positionProbabilityMatrix

    /// Increase counter of PositionProbabilityMatrix at fixed position by n.
    let increaseInPlacePPMBy (pos:int) (bioItem:'a when 'a :> IBioItem) n (positionProbabilityMatrix:PositionProbabilityMatrix) = 
        positionProbabilityMatrix.[bioItem, pos] <- positionProbabilityMatrix.[bioItem, pos] + n
        positionProbabilityMatrix

    /// Create new PositionWeightMatrix based on existing PositionFrequencyMatrix. 
    /// The counts of each position of each element are transformed to floats.
    let createPPMOf (positionFrequencyMatrix:PositionFrequencyMatrix) =
        positionFrequencyMatrix.Matrix |> Array2D.map (fun item -> float item)
        |> fun item -> new PositionProbabilityMatrix(item)

    // Create new PositionWeightMatrix based on existing PositionFrequencyMatrix. 
    /// The counts of each position of each element are transformed to floats.
    let normalizePPM (sourceCount:int) (alphabet:#IBioItem[]) (pseudoCount:float) (positionFrequencyMatrix:PositionProbabilityMatrix) =
        let positionProbabilityMatrix = new PositionProbabilityMatrix(positionFrequencyMatrix.Matrix)
        let sum = (float sourceCount) + ((float alphabet.Length) * pseudoCount)
        for item in alphabet do
            for position = 0 to (Array2D.length2 positionProbabilityMatrix.Matrix) - 1 do
            positionProbabilityMatrix.[item, position] <- (positionProbabilityMatrix.[item, position] + pseudoCount)/sum
        positionProbabilityMatrix
    
    type PositionWeightMatrix internal(matrix:float [,]) =

        inherit BaseMatrix<IBioItem, float>(matrix)

        new(rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionWeightMatrix(arr)

    /// Increase counter of PositionWeightMatrix at fixed position by 1.
    let increaseInPlacePWM (pos:int) (bioItem:'a when 'a :> IBioItem) (positionWeightMatrix:PositionWeightMatrix) = 
        positionWeightMatrix.[bioItem, pos] <- positionWeightMatrix.[bioItem, pos] + 1.
        positionWeightMatrix

    // Increase counter of PositionWeightMatrix at fixed position by n.
    let increaseInPlacePWMBy (bioItem:'a when 'a :> IBioItem) (pos:int) n (positionWeightMatrix:PositionWeightMatrix) = 
        positionWeightMatrix.[bioItem, pos] <- positionWeightMatrix.[bioItem, pos] + n
        positionWeightMatrix

    /// Create PositionWeightMatrix based on ProbabilityCompositeVector and PositionProbabilityMatrix.
    let createPositionWeightMatrix (alphabet:#IBioItem[]) (pcv:ProbabilityCompositeVector) (ppMatrix:PositionProbabilityMatrix) =
        let pwMatrix = new PositionWeightMatrix(Array2D.length2 ppMatrix.Matrix)
        for item in alphabet do
            for position=0 to (Array2D.length2 ppMatrix.Matrix)-1 do
                pwMatrix.[item, position] <- ppMatrix.[item, position]/pcv.[item]
        pwMatrix

    /// Calculate the score of the given sequence based on the Probabilities of the PositionWeightMatrix.
    let calculateSegmentScoreBy (pwMatrix:PositionWeightMatrix) (bioItems:BioArray.BioArray<#IBioItem>) =
        Array.fold (fun (position:int, value:float) (bios:#IBioItem) -> 
            position + 1, value * (pwMatrix.[bios, position])) (0, 1.) bioItems
        |> snd
