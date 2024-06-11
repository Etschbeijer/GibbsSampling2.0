module HelperFunctions

open FSharpAux

/// Replace element at index
let replaceIndexElement (result:array<'a>) (targetI:array<'a>) (targetII:array<'a>) (index:int) = 
    
    if index < targetI.Length then
        result.[index] <- targetI.[index]
    if index < targetII.Length then
        result.[targetI.Length + index] <- targetII[index]
    result
    
/// Merge Arrays
let mergeArrays (targetI:array<'a>) (targetII:array<'a>) =
    
    if targetI.Length > targetII.Length then
        targetI
        |> Array.foldi (fun i (element) acc -> 
            (replaceIndexElement element targetI targetII i))
                (Array.zeroCreate (targetI.Length + targetII.Length))
    else
        targetII
        |> Array.foldi (fun i (element) acc -> 
            (replaceIndexElement element targetI targetII i))
                (Array.zeroCreate (targetI.Length + targetII.Length))
