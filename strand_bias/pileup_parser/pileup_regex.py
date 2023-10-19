
import re

READ_START = re.compile(r"\^.")  # Start symbol (^) followed by read mapping quality.

MATCH = r"\.,"  # Base matches the reference base.
MISMATCH = "ACGTNacgtn"  # Base is a mismatch to the reference base.
SKIP = "><"  # Reference skip (due to CIGAR "N").
DELETION = r"\*#"  # Deletion of the reference base (CIGAR "D").

INDEL = re.compile(rf"[\+-]\d+[ACGTNacgtn\*#]+")
READ_END = re.compile("\$")

PILEUP_POS_RE = re.compile(rf"({READ_START.pattern})?"
                           rf"([{MATCH}{MISMATCH}{SKIP}{DELETION}])"
                           rf"({INDEL.pattern})?"
                           rf"({READ_END.pattern})?")

FORWARD_BASES = set(".ACGTN>*")
REVERSE_BASES = set(",acgtn<#")
MATCH_SET = set(".,")
