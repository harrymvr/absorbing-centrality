class CanonicalizationError(Exception):
    """Exception related to the graph canonicalization procedure.
    """

    def __init__(self, message):
        """Initialize the exception and supply a descriptive message.

        Parameters
        ----------
        message : str
            The description of the error.
        """

        super(CanonicalizationError, self). __init__(message)
