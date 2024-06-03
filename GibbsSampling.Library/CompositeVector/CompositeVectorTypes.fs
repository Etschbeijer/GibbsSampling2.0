namespace CompositeVector

open System
open BioFSharp

module Types =
    
    /// Matrix with fixed positions for nucleotides and amino acids.
    type BaseMatrix<'a, 'value when 'a :> IBioItem>internal(matrix:'value [,]) =

        let getRowArray2DIndex (key:'a) =
            (BioItem.symbol key |> int) - 42

        new (rowLength:int) =
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new BaseMatrix<_, 'value>(arr)

        member val _matrix = 
            if (obj.ReferenceEquals(matrix, null)) then
                raise (ArgumentNullException("array"))        
            (Array2D.copy matrix) with get, set

        member this.Matrix = this._matrix

        member this.Item
            with get (column, row)       = this._matrix.[getRowArray2DIndex column, row]
            and  set (column, row) value = this._matrix.[getRowArray2DIndex column, row] <- value
    
    /// One dimensional array with fixed positions for each element.
    type CompositeVector<'a, 'value when 'a :> IBioItem>internal (array:'value []) =
    
        let getIndex (key:'a) =
            (BioItem.symbol key |> int) - 42

        new () =
            let arr:'value [] = Array.zeroCreate 49
            new CompositeVector<_, 'value>(arr)

        member this.Array = 
            if (obj.ReferenceEquals(array, null)) then
                raise (ArgumentNullException("array"))
            array
            
        member this.Item
            with get i       = array.[getIndex i]
            and  set i value = array.[getIndex i] <- value

    /// One dimensional array with fixed positions for each element.
    /// Use to track frequency of elements independent of position in source.
    type FrequencyCompositeVector internal (array:int []) =

        inherit CompositeVector<IBioItem, int>(array)

        new() = 
            let arr:'value [] = Array.zeroCreate 49
            new FrequencyCompositeVector(arr)

    /// Matrix with fixed positions for nucleotides and amino acids with default value of 0 for probability.
    type PositionProbabilityMatrix internal(matrix:float [,]) =

        inherit BaseMatrix<IBioItem, float>(matrix)

        new(rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionProbabilityMatrix(arr)

    /// One dimensional array with fixed positions for each element.
    /// Use to track probability of elements independent of position in source.
    type ProbabilityCompositeVector internal (array:float []) =
        
        inherit CompositeVector<IBioItem,float>(array)

        new() = 
            let arr:'value [] = Array.zeroCreate 49
            new ProbabilityCompositeVector(arr)
