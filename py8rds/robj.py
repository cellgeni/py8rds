import numpy as np

INT_NA = -2147483648


class RdsFile:
    def __init__(self):
        self.format_version = 0
        self.writer_version = [0, 0, 0]
        self.reader_version = [0, 0, 0]
        self.encoding = ""
        self.object = None
        self.environments = []
        self.symbols = []
        self.external_pointers = []


class Robj:
    def __init__(self):
        # can be: Robj or list of scalars or Robjs or list of [scalar,Robj] pairs/Nones
        self.value = None
        # attributes has unique names and then can be stored in dict, but they are serialized as pairlist which names are not necessary unique.
        # So we will store attributes as list of [scalar,Robj] pairs
        self.attributes = []

    def get(self, inxs):
        """
        Recursively navigate into an Robj using indices and/or names.
        Indices are used only to search in the value; names are used to search both in attributes (first)
        and then in values if they are a pairlist (should be the case for S4 slots).

        Parameters
        ----------
        inxs : int | str | list[int | str]
            A single index/name or a list describing a path within an R object.

        Returns
        -------
        Any
            - Another `Robj`,
            - a primitive Python value (int, float, bool, str, etc.),
            - or `None` if any step along the path cannot be resolved.

        Examples
        --------
        Assuming an R list like `list(a = 1:3, b = 4:5)`:

        >>> robj.get("a")
        <Robj for 1:3>
        >>> robj.get(["a", 1])
        2
        """
        if not isinstance(inxs, list):
            inxs = [inxs]
        r = self._get(inxs[0])
        if (r is None) or (len(inxs) == 1):
            return r
        return r.get(inxs[1:])

    def _get(self, inx):
        # look by index
        if isinstance(inx, int):
            if isinstance(self.value, list):
                return self.value[inx]
            elif inx == 0:
                return self.value
            else:
                return None
        # look by name
        elif isinstance(inx, str):
            # in attributes
            for a in self.attributes:
                if inx == a[0]:
                    return a[1]
            # in values
            if not isinstance(self.value, list):
                return None
            for v in self.value:
                if isinstance(v, list):
                    if inx == v[0]:
                        return v[1]
        return None

    def getClass(self):
        r = self._get("class")
        if r is None:
            r = self._get("sexptype")
        else:
            r = ",".join(r.value)
        return r

    def is_primitive(self, x):
        return isinstance(x, (int, float, bool, str, bytes, type(None), np.generic))

    def show(self, level=1):
        print(self.toString(level=level))

    def toString(self, name="Robj", indent="", maxItems=5, level=1e3):
        if level < 0:
            return ""
        level -= 1
        parts = []
        parts.append(indent + name + "(" + str(self.getClass()) + "): ")
        indent = indent.replace("+", "|").replace("*", "|").replace("&", "|")
        # value
        values = self.value
        values_is_array = isinstance(values, np.ndarray)
        if (not isinstance(values, list)) and (not values_is_array):
            values = [values]
        handled_values = False
        if values_is_array:
            if values.size == 0:
                parts.append("[]\n")
                handled_values = True
            else:
                first_val = values.flat[0]
                if self.is_primitive(first_val):
                    parts.append("[")
                    preview_count = min(3, values.size)
                    parts.append(
                        ",".join([str(values.flat[i]) for i in range(preview_count)])
                    )
                    if values.size > 3:
                        parts.append(",...")
                    parts.append("]")
                    if values.ndim > 1:
                        parts.append(f" shape={values.shape}")
                    parts.append("\n")
                    handled_values = True
                else:
                    values = values.tolist()
                    values_is_array = False
        if not handled_values:
            if len(values) == 0:
                parts.append("[]\n")
            elif self.is_primitive(values[0]):
                parts.append("[")
                parts.append(",".join([str(v) for v in values[: min(3, len(values))]]))
                if len(values) > 3:
                    parts.append(",...")
                parts.append("]\n")
            elif isinstance(values[0], Robj):
                parts.append("\n")
                for i in range(len(values)):
                    parts.append(
                        values[i].toString(
                            name=str(i),
                            indent=indent + "+",
                            maxItems=maxItems,
                            level=level,
                        )
                    )
            else:
                parts.append("\n")
                for i in range(len(values)):
                    if self.is_primitive(values[i]):
                        parts.append(indent + "&" + str(values[i]) + "\n")
                    elif isinstance(values[i][1], Robj):
                        parts.append(
                            values[i][1].toString(
                                name=str(values[i][0]),
                                indent=indent + "&",
                                maxItems=maxItems,
                                level=level,
                            )
                        )
                    else:
                        parts.append(
                            indent
                            + "&"
                            + str(values[i][0])
                            + ":"
                            + str(values[i][1])
                            + "\n"
                        )
        # attributes
        for a in self.attributes:
            if len(a) != 2:
                raise RuntimeError("Attributes are not in a pair: " + str(a))
            if a[0] != "sexptype":
                parts.append(
                    a[1].toString(
                        name=str(a[0]),
                        indent=indent + "*",
                        maxItems=maxItems,
                        level=level,
                    )
                )
        return "".join(parts)
