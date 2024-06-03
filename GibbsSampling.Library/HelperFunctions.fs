module HelperFunctions

open FSharpAux

/// Replace element at index
let replaceIndexElement (result:array<'a>) (targetsILength:int) (targetI:'a)(targetsII:array<'a>) (index:int) = 
        result.[index] <- targetI
        if index < targetsII.Length then
            result.[targetsILength + index] <- targetsII[index]
        result
    
/// Merge Arrays
let mergeArrays (targetI:array<'a>) (targetII:array<'a>) =
    if targetI.Length > targetII.Length then
        targetI
        |> Array.foldi (fun i (element) acc -> 
            (replaceIndexElement element targetI.Length acc targetII i))
                (Array.zeroCreate (targetI.Length + targetII.Length))
    else
        targetII
        |> Array.foldi (fun i (element) acc -> 
            (replaceIndexElement element targetII.Length acc targetI i))
                (Array.zeroCreate (targetI.Length + targetII.Length))