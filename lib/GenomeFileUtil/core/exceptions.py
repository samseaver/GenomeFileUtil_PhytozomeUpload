import json


class RENotFound(RuntimeError):
    """A resource we need from the Relation Engine was not found."""

    def __init__(self, resource, key, val, resp_json):
        """
        `key` - the key we used to try to find something.
        `val` - the val we used to try to look up the above key.
        """
        self.resource = resource
        self.key = key
        self.val = val
        self.resp_json = resp_json

    def __repr__(self):
        # Eg: "The Relation Engine was unable to fetch taxon results by taxonomy ID
        #     using the value 123. The server response was:
        #     etc."
        return (f"The Relation Engine was unable to fetch {self.resource} results by "
                f"{self.key} using the value '{self.val}'. The server response "
                f"was: \n{json.dumps(self.resp_json, indent=2)}")
