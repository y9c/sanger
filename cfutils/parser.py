#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Copyright 2011 by Wibowo Arindrarto (w.arindrarto@gmail.com)
# Revisions copyright 2011 by Peter Cock.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
#
# Created: 2019-08-30 13:51


"""parser for the ABI format.

Learned from Bio.SeqIO 

ABI is the format used by Applied Biosystem's sequencing machines to store
sequencing results.

For more details on the format specification, visit:
http://www.appliedbiosystem.com/support/software_community/ABIF_File_Format.pdf
"""

import datetime
import struct
from pathlib import Path

# dictionary for determining which tags goes into SeqRecord annotation
# each key is tag_name + tag_number
# if a tag entry needs to be added, just add its key and its key
# for the annotations dictionary as the value
_EXTRACT = {
    "TUBE1": "sample_well",
    "DySN1": "dye",
    "GTyp1": "polymer",
    "MODL1": "machine_model",
}
# dictionary for tags that require preprocessing before use in creating
# seqrecords
_SPCTAGS = [
    "PBAS2",  # base-called sequence
    "PCON2",  # quality values of base-called sequence
    "SMPL1",  # sample id inputted before sequencing run
    "RUND1",  # run start date
    "RUND2",  # run finish date
    "RUNT1",  # run start time
    "RUNT2",  # run finish time
    # NOTE: The following are used for trace data
    "PLOC2",  # position of peaks
    "DATA1",  # channel1 raw data
    "DATA2",  # channel2 raw data
    "DATA3",  # channel3 raw data
    "DATA4",  # channel4 raw data
    "DATA9",  # channel1 analyzed data
    "DATA10",  # channel2 analyzed data
    "DATA11",  # channel3 analyzed data
    "DATA12",  # channel4 analyzed data
    "FWO_1",  # base order for channels
]
# dictionary for data unpacking format
_BYTEFMT = {
    1: "b",  # byte
    2: "s",  # char
    3: "H",  # word
    4: "h",  # short
    5: "i",  # long
    6: "2i",  # rational, legacy unsupported
    7: "f",  # float
    8: "d",  # double
    10: "h2B",  # date
    11: "4B",  # time
    12: "2i2b",  # thumb
    13: "B",  # bool
    14: "2h",  # point, legacy unsupported
    15: "4h",  # rect, legacy unsupported
    16: "2i",  # vPoint, legacy unsupported
    17: "4i",  # vRect, legacy unsupported
    18: "s",  # pString
    19: "s",  # cString
    20: "2i",  # tag, legacy unsupported
}
# header data structure (exluding 4 byte ABIF marker)
_HEADFMT = ">H4sI2H3I"
# directory data structure
_DIRFMT = ">4sI2H4I"


class SeqRecord:
    def __init__(self, seq, id="", name="", description="", annotations=None, letter_annotations=None):
        self.seq = str(seq)  # Ensure sequence is stored as a string
        self.id = id
        self.name = name
        self.description = description
        self.annotations = annotations if annotations is not None else {}
        self.letter_annotations = letter_annotations if letter_annotations is not None else {}

    def __getitem__(self, key):
        new_seq = self.seq[key]
        new_annotations = self.annotations.copy()
        new_letter_annotations = {k: v[key] for k, v in self.letter_annotations.items()}
        return SeqRecord(new_seq, self.id, self.name, self.description, new_annotations, new_letter_annotations)

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq  # Return the sequence string for display


def abi_iterator(handle):
    """Iterator for the Abi file format."""

    # raise exception if handle mode is not 'rb'
    if hasattr(handle, "mode"):
        if set("rb") != set(handle.mode.lower()):
            raise ValueError("ABI files has to be opened in 'rb' mode.")

    # check if input file is a valid Abi file
    handle.seek(0)
    marker = handle.read(4)
    if not marker:
        # handle empty file gracefully
        raise StopIteration
    if marker != b"ABIF":
        raise IOError("File should start ABIF, not %r" % marker)

    # dirty hack for handling time information
    times = {"RUND1": "", "RUND2": "", "RUNT1": "", "RUNT2": ""}

    # initialize annotations
    annot = dict(zip(_EXTRACT.values(), [None] * len(_EXTRACT)))

    # parse header and extract data from directories
    header = struct.unpack(_HEADFMT, handle.read(struct.calcsize(_HEADFMT)))

    for tag_name, tag_number, tag_data in _abi_parse_header(header, handle):
        # stop iteration if all desired tags have been extracted
        # 4 tags from _EXTRACT + 2 time tags from _SPCTAGS - 3,
        # and seq, qual, id
        # todo

        key = tag_name + str(tag_number)

        # PBAS2 is base-called sequence
        if key == "PBAS2":
            seq = tag_data
        # PCON2 is quality values of base-called sequence
        elif key == "PCON2":
            qual = [ord(val) for val in tag_data]
        # PLOC2 is the location of peaks
        elif key == "PLOC2":
            peakamps = [float(val) for val in tag_data]
            annot["peak positions"] = peakamps
        # DATA1-DATA4 is raw channel 1-4 output, DATA9-12 the analyzed one
        elif key in ["DATA9", "DATA10", "DATA11", "DATA12"]:
            rawch = [float(val) for val in tag_data]
            annot["channel " + str(int(key[4:]) - 8)] = rawch
        # FWO_1 is the order of channels in bases
        elif key == "FWO_1":
            channelorders = tag_data
            annot["channels"] = channelorders
        # SMPL1 is sample id entered before sequencing run
        elif key == "SMPL1":
            sample_id = tag_data
        elif key in times:
            times[key] = tag_data
        else:
            # extract sequence annotation as defined in _EXTRACT
            if key in _EXTRACT:
                annot[_EXTRACT[key]] = tag_data

    # set time annotations
    annot["run_start"] = "%s %s" % (times["RUND1"], times["RUNT1"])
    annot["run_finish"] = "%s %s" % (times["RUND2"], times["RUNT2"])

    # use the file name as SeqRecord.name if available
    #  try:
    file_name = Path(handle.name).stem
    #  except:
    #  file_name = ""

    record = SeqRecord(
        seq,  # Use the sequence string directly
        id=sample_id,
        name=file_name,
        description="",
        annotations=annot,
        letter_annotations={"phred_quality": qual},
    )

    yield record


def _abi_parse_header(header, handle):
    """Generator that returns directory contents."""
    # header structure (after ABIF marker):
    # file version, tag name, tag number,
    # element type code, element size, number of elements
    # data size, data offset, handle (not file handle)
    head_elem_size = header[4]
    head_elem_num = header[5]
    head_offset = header[7]
    index = 0

    while index < head_elem_num:
        start = head_offset + index * head_elem_size
        # add directory offset to tuple
        # to handle directories with data size <= 4 bytes
        handle.seek(start)
        dir_entry = struct.unpack(_DIRFMT, handle.read(struct.calcsize(_DIRFMT))) + (
            start,
        )
        index += 1
        # only parse desired dirs
        key = dir_entry[0].decode()
        key += str(dir_entry[1])
        if key in list(_EXTRACT.keys()) + _SPCTAGS:
            tag_name = dir_entry[0].decode()
            tag_number = dir_entry[1]
            elem_code = dir_entry[2]
            elem_num = dir_entry[4]
            data_size = dir_entry[5]
            data_offset = dir_entry[6]
            tag_offset = dir_entry[8]
            # if data size <= 4 bytes, data is stored inside tag
            # so offset needs to be changed
            if data_size <= 4:
                data_offset = tag_offset + 20
            handle.seek(data_offset)
            data = handle.read(data_size)
            yield tag_name, tag_number, _parse_tag_data(elem_code, elem_num, data)


def _abi_trim(seq_record: SeqRecord) -> SeqRecord:
    """Trims the sequence using Richard Mott's modified trimming algorithm.

    seq_record - SeqRecord object to be trimmed.

    Trimmed bases are determined from their segment score, which is a
    cumulative sum of each base's score. Base scores are calculated from
    their quality values.

    More about the trimming algorithm:
    http://www.phrap.org/phredphrap/phred.html
    http://www.clcbio.com/manual/genomics/Quality_abif_trimming.html
    """

    start = False  # flag for starting position of trimmed sequence
    segment = 20  # minimum sequence length
    trim_start = 0  # init start index
    cutoff = 0.05  # default cutoff value for calculating base score

    if len(seq_record) <= segment:
        return seq_record
    else:
        # calculate base score
        score_list = [
            cutoff - (10 ** (qual / -10.0))
            for qual in seq_record.letter_annotations["phred_quality"]
        ]

        # calculate cummulative score
        # if cummulative value < 0, set it to 0
        # first value is set to 0, because of the assumption that
        # the first base will always be trimmed out
        cummul_score = [0]
        for i in range(1, len(score_list)):
            score = cummul_score[-1] + score_list[i]
            if score < 0:
                cummul_score.append(0)
            else:
                cummul_score.append(score)
                if not start:
                    # trim_start = value when cummulative score is first > 0
                    trim_start = i
                    start = True

        # trim_finish = index of highest cummulative score,
        # marking the end of sequence segment with highest cummulative score
        trim_finish = cummul_score.index(max(cummul_score))
        new_record = seq_record[trim_start:trim_finish]
        new_record.annotations = seq_record.annotations.copy()

        return new_record


def _parse_tag_data(elem_code, elem_num, raw_data):
    """Returns single data value.

    elem_code - What kind of data
    elem_num - How many data points
    raw_data - abi file object from which the tags would be unpacked
    """
    if elem_code in _BYTEFMT:
        # because '>1s' unpack differently from '>s'
        if elem_num == 1:
            num = ""
        else:
            num = str(elem_num)
        fmt = ">" + num + _BYTEFMT[elem_code]

        assert len(raw_data) == struct.calcsize(fmt)
        data = struct.unpack(fmt, raw_data)

        # no need to use tuple if len(data) == 1
        # also if data is date / time
        if elem_code not in [10, 11] and len(data) == 1:
            data = data[0]

        # account for different data types
        if elem_code == 2:
            return data.decode()
        if elem_code == 10:
            return str(datetime.date(*data))
        if elem_code == 11:
            return str(datetime.time(*data[:3]))
        if elem_code == 13:
            return bool(data)
        if elem_code == 18:
            return data[1:].decode()
        if elem_code == 19:
            return data[:-1].decode()
        return data
    else:
        return None


def trim_and_rescale_trace(seq):
    """Trim traces to peak positions, shift to start from zero, and rescale."""

    traces = [seq.annotations["channel " + str(i)] for i in range(1, 5)]
    peaks = seq.annotations["peak positions"]
    n = len(peaks)
    step = 1.0 * (peaks[-1] - peaks[0]) / n

    traces = [
        [t for (i, t) in enumerate(trace) if peaks[0] <= i < peaks[-1]]
        for trace in traces
    ]
    peaks = [(p - peaks[0]) / step for p in peaks]

    x = [1.0 * i / step for i in range(len(traces[0]))]

    seq.annotations["peak positions"] = peaks
    for i, trace in enumerate(traces, 1):
        seq.annotations["channel " + str(i)] = trace
    seq.annotations["trace_x"] = x
    return seq


def rescale_trace(seq: SeqRecord) -> SeqRecord:
    traces = [seq.annotations["channel " + str(i)] for i in range(1, 5)]
    peaks = seq.annotations["peak positions"]
    n = len(peaks)
    step = 1.0 * (peaks[-1] - peaks[0]) / n
    traces = [
        [t for (i, t) in enumerate(trace) if peaks[0] <= i < peaks[-1]]
        for trace in traces
    ]
    #  peaks = [(p - peaks[0]) / step for p in peaks]
    peaks = [p / step for p in peaks]

    x = [1.0 * i / step for i in range(len(traces[0]))]

    seq.annotations["peak positions"] = peaks
    for i, trace in enumerate(traces, 1):
        seq.annotations["channel " + str(i)] = trace
    seq.annotations["trace_x"] = x
    return seq


def parse_abi(filename: str, rescale: bool = True) -> SeqRecord:
    """Parse an ABI file from Sanger sequencing.

    Args:
        filename: path to an ``.ab1`` file.
        rescale: if True (default), rescale the trace coordinates so the
            peak axis is in 0-based sequence positions (convenient for
            plotting).  If False, keep the raw trace arrays and sample-index
            peak positions, which is required for peak/baseline analysis such
            as :func:`cfutils.basecaller.call_bases`.
    """
    with open(filename, "rb") as abifile:
        seq = list(abi_iterator(abifile))[0]

    if rescale:
        seq = rescale_trace(seq)
    return seq


def parse_fasta(filename: str) -> SeqRecord:
    """parse_fasta may support other type of file in the future.

    :param filename:
    :type filename: str
    :rtype: SeqRecord
    """
    with open(filename, "r") as file:
        lines = file.readlines()
        id_line = lines[0].strip()
        sequence = ''.join(line.strip() for line in lines[1:])
        seq_id = id_line[1:] if id_line.startswith('>') else ''
        return SeqRecord(sequence, id=seq_id)
