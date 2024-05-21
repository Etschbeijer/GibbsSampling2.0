namespace MotifSampler

module Types =

    /// Contains the information of the probability and the start position(s) of the found motifs(s). 
    type MotifIndex =
        {
            PWMS        : float
            Positions   : int list
        }