namespace PositionMatrix

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
        
        member this.Matrix = 
            if (obj.ReferenceEquals(matrix, null)) then
                raise (ArgumentNullException("array"))
            matrix
            
        member this.Item
            with get (column, row)       = matrix.[getRowArray2DIndex column, row]
            and  set (column, row) value = matrix.[getRowArray2DIndex column, row] <- value

    /// Matrix with fixed positions for nucleotides and amino acids with default value of 0 for frequency.
    type PositionFrequencyMatrix internal(matrix:int [,]) =

        inherit BaseMatrix<IBioItem, int>(matrix)

        new (rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionFrequencyMatrix(arr)

    /// Matrix with fixed positions for nucleotides and amino acids with default value of 0 for probability.
    type PositionProbabilityMatrix internal(matrix:float [,]) =

        inherit BaseMatrix<IBioItem, float>(matrix)

        new(rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionProbabilityMatrix(arr)

    type PositionWeightMatrix internal(matrix:float [,]) =

        inherit BaseMatrix<IBioItem, float>(matrix)

        new(rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionWeightMatrix(arr)
