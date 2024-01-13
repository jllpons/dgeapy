from dgeapy_utils import Color, TermMsg


class InvalidArgumentError(Exception):
    """Raised when an invalid argument is passed to the program."""

    def __init__(self, invalid_argument,
                 message=f"{TermMsg.ERROR}: Invalid argument passed to the program",
                 closure=""):
        self.invalid_argument = invalid_argument
        self.message = message
        self.closure = closure

        super().__init__(self.message)

    def __str__(self):
        return (f"{self.message}: {Color.YELLOW}{self.invalid_argument}{Color.RESET}"
                + self.closure)


class FileNotFoundError(Exception):
    """Raised when a file is not found."""

    def __init__(self, file_name, message=f"{TermMsg.ERROR}: File not found"):
        self.file_name = file_name
        self.message = message

        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}: {Color.YELLOW}{self.file_name}{Color.RESET}"


class MissingArgumentError(Exception):
    """Raised when a required argument is missing."""

    def __init__(self, missing_argument,
                 message=f"{TermMsg.ERROR}: Missing required argument", closure=""):
        self.missing_argument = missing_argument
        self.message = message
        self.closure = closure

        super().__init__(self.message)

    def __str__(self):
        return (f"{self.message}: {Color.YELLOW}{self.missing_argument}{Color.RESET}"
                f"{self.closure}")


class UnsupportedFileFormatError(Exception):
    """Raised when a file format is not supported."""

    def __init__(self, file_format,
                 message=f"{TermMsg.ERROR}: File format not supported", closure=""):
        self.file_format = file_format
        self.message = message
        self.closure = closure

        super().__init__(self.message)

    def __str__(self):
        return (f"{self.message}: {Color.YELLOW}{self.file_format}{Color.RESET}"
                f"{self.closure}")


class RequiredColumnNotFoundError(Exception):
    """Raised when a required column is not found in a table."""

    def __init__(self, column_name, column_type, message=None, closure=""):
        self.column_name = column_name
        self.column_type = column_type
        self.closure = closure

        # Set the default message if none is provided
        if message is None:
            message = f"{TermMsg.ERROR}: Required {Color.YELLOW} {self.column_type} {Color.RESET} column not found"

        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return (f"{self.message}, given name was: {Color.YELLOW}{self.column_name}{Color.RESET}"
                f"{self.closure}")


class NanValuesInIndexColumnError(Exception):
    """Raised when there are NaN values in the index column."""

    def __init__(self, message=f"{TermMsg.ERROR}: NaN values in index column",
                 closure=""):
        self.message = message
        self.closure = closure

        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}{self.closure}"


