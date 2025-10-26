class InputParser:
    """Parser for custom input files with 'arg : value' format."""

    def __init__(self, filename=None):
        """
        Initialize the parser.

        Args:
            filename (str, optional): Path to the input file to parse immediately
        """
        self.data = {}
        if filename:
            self.parse_file(filename)

    def parse_file(self, filename):
        """
        Parse an input file and store the arguments.

        Args:
            filename (str): Path to the input file

        Returns:
            dict: Dictionary of parsed arguments
        """
        self.data = {}
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                # Skip comment lines
                if line.startswith('#'):
                    continue

                # Remove inline comments (everything after #)
                if '#' in line:
                    line = line.split('#')[0].strip()

                # Skip if line became empty after removing comment
                if not line:
                    continue

                # Check if line contains ':'
                if ':' not in line:
                    print(f"Warning: Line {line_num} has no ':' separator, skipping")
                    continue

                # Split by ':' and clean up whitespace
                parts = line.split(':', 1)  # Split only on first ':'
                arg = parts[0].strip()
                value = parts[1].strip()

                # Parse and convert the value
                self.data[arg] = self._parse_value(value)

        return self.data

    def _parse_value(self, value_str):
        """
        Parse and convert a value string to appropriate type.
        Handles lists (separated by ,) and automatic type conversion.

        Args:
            value_str (str): The value string to parse

        Returns:
            Converted value (int, float, str, or list of these types)
        """
        # Check if it's a list (contains semicolon)
        if ',' in value_str:
            # Split by semicolon and parse each value
            values = [v.strip() for v in value_str.split(',')]
            return [self._convert_type(v) for v in values if v]  # Filter empty strings
        else:
            # Single value - convert type
            return self._convert_type(value_str)

    def _convert_type(self, value_str):
        """
        Automatically convert a string to its appropriate type.

        Args:
            value_str (str): String to convert

        Returns:
            int, float, or str depending on the content
        """
        # Try integer first
        try:
            # Check if it looks like an integer (no decimal point)
            if '.' not in value_str:
                return int(value_str)
        except ValueError:
            pass

        # Try float
        try:
            return float(value_str)
        except ValueError:
            pass

        # Return as string if nothing else works
        return value_str

    def parse_string(self, text):
        """
        Parse a string with multiple lines.

        Args:
            text (str): Multi-line string to parse

        Returns:
            dict: Dictionary of parsed arguments
        """
        self.data = {}
        for line_num, line in enumerate(text.split('\n'), 1):
            line = line.strip()

            if not line:
                continue

            # Skip comment lines
            if line.startswith('#'):
                continue

            # Remove inline comments (everything after #)
            if '#' in line:
                line = line.split('#')[0].strip()

            # Skip if line became empty after removing comment
            if not line:
                continue

            if ':' not in line:
                print(f"Warning: Line {line_num} has no ':' separator, skipping")
                continue

            parts = line.split(':', 1)
            arg = parts[0].strip()
            value = parts[1].strip()

            # Parse and convert the value
            self.data[arg] = self._parse_value(value)

        return self.data

    def get(self, arg, default=None):
        """
        Get the value for a specific argument.

        Args:
            arg (str): Argument name
            default: Default value if argument not found

        Returns:
            The value associated with the argument, or default if not found
        """
        return self.data.get(arg, default)

    def get_int(self, arg, default=None):
        """Get argument value as integer."""
        value = self.get(arg)
        if value is None:
            return default
        try:
            return int(value)
        except ValueError:
            print(f"Warning: Cannot convert '{value}' to int for arg '{arg}'")
            return default

    def get_float(self, arg, default=None):
        """Get argument value as float."""
        value = self.get(arg)
        if value is None:
            return default
        try:
            return float(value)
        except ValueError:
            print(f"Warning: Cannot convert '{value}' to float for arg '{arg}'")
            return default

    def get_bool(self, arg, default=None):
        """Get argument value as boolean (true/false, yes/no, 1/0)."""
        value = self.get(arg)
        if value is None:
            return default

        value_lower = value.lower()
        if value_lower in ('true', 'yes', '1', 'on'):
            return True
        elif value_lower in ('false', 'no', '0', 'off'):
            return False
        else:
            print(f"Warning: Cannot convert '{value}' to bool for arg '{arg}'")
            return default

    def __getitem__(self, arg):
        """Allow dictionary-style access: parser['arg']"""
        return self.data[arg]

    def __contains__(self, arg):
        """Allow 'in' operator: 'arg' in parser"""
        return arg in self.data

    def __str__(self):
        """String representation of parsed data."""
        return str(self.data)

    def items(self):
        """Return items like a dictionary."""
        return self.data.items()

    def keys(self):
        """Return keys like a dictionary."""
        return self.data.keys()

    def values(self):
        """Return values like a dictionary."""
        return self.data.values()


# Example usage
if __name__ == "__main__":
    # Example 1: Parse from string with automatic type conversion and lists
    sample_input = """
    # This is a configuration file
    name : John Doe  # User's full name (string)
    age    :    30   # Integer
    email:john@example.com  # Contact email

    # Status settings
    active : true
    score   :   95.5  # Float

    # List examples (separated by semicolons)
    coordinates : 10, 20, 30
    prices : 19.99, 29.99, 39.99
    colors : red, green, blue
    mixed : 100, 3.14, text, 42
    """

    parser = InputParser()
    parser.parse_string(sample_input)

    print("Parsed data:")
    for key, value in parser.items():
        print(f"{key}: {value} (type: {type(value).__name__})")
    print()

    # Access values
    print(f"Name: {parser.get('name')} (type: {type(parser.get('name')).__name__})")
    print(f"Age: {parser.get('age')} (type: {type(parser.get('age')).__name__})")
    print(f"Score: {parser.get('score')} (type: {type(parser.get('score')).__name__})")
    print()

    # Access list values
    print(f"Coordinates: {parser.get('coordinates')}")
    print(f"Prices: {parser.get('prices')}")
    print(f"Colors: {parser.get('colors')}")
    print(f"Mixed list: {parser.get('mixed')}")
    print()

    # Example 2: Parse from file
    # parser = InputParser('config.txt')
    # print(parser.get('some_arg'))
