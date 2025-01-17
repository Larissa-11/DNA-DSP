import copy

from Chamaeleo.utils import data_handle, indexer
from Chamaeleo.utils.pipelines import TranscodePipeline

from utils.my_chamaeleo import data_handle as my_data_handle


class MyTranscodePipeline(TranscodePipeline):
    def transcode(
            self,
            direction=None,
            input_path=None,
            input_byte=None,
            input_string=None,
            index_length=None,
            segment_length=120,
            **info
    ):
        if direction is not None:
            if direction == "t_c":
                self.records["payload length"] = segment_length

                if input_path is not None:
                    bit_segments, bit_size = data_handle.read_bits_from_file(input_path, segment_length,
                                                                             self.need_logs)
                elif input_byte is not None:
                    bit_segments, bit_size = my_data_handle.read_bits_from_byte(input_byte, segment_length,
                                                                                self.need_logs)
                elif input_string is not None:
                    bit_segments, bit_size = data_handle.read_bits_from_str(input_string, segment_length,
                                                                            self.need_logs)
                else:
                    raise ValueError("There is no digital data input here!")

                if info.get("index"):
                    if index_length:
                        bit_segments, index_length = indexer.connect_all(bit_segments, index_length, self.need_logs)
                    else:
                        bit_segments, index_length = indexer.connect_all(bit_segments, None, self.need_logs)

                    self.records["index length"] = index_length
                else:
                    self.records["index length"] = 0

                if self.error_correction is not None:
                    bit_segments, error_correction_length = self.error_correction.insert(bit_segments)
                    self.records["error-correction length"] = error_correction_length
                else:
                    self.records["error-correction length"] = 0

                results = self.coding_scheme.silicon_to_carbon(bit_segments, bit_size)

                dna_sequences = results["dna"]

                self.records["information density"] = round(results["i"], 3)
                self.records["encoding runtime"] = round(results["t"], 3)

                if "output_path" in info:
                    data_handle.write_dna_file(info["output_path"], dna_sequences, self.need_logs)

                return {
                    "bit": bit_segments,
                    "dna": dna_sequences,
                    "info": {
                        "bit_size": self.coding_scheme.bit_size,
                        "segment_length": self.coding_scheme.segment_length,
                        "total_count": getattr(self.coding_scheme, "total_count",
                                               getattr(self.coding_scheme, "decode_packets", 0)),
                        "index_length": getattr(self.coding_scheme, "index_length", 0)
                    }
                }
            elif direction == "t_s":
                if input_path is not None:
                    dna_sequences = data_handle.read_dna_file(input_path, self.need_logs)
                elif input_byte is not None:
                    dna_sequences = my_data_handle.read_dna_byte(input_byte, self.need_logs)
                elif input_string is not None:
                    dna_sequences = []
                    for index, string in enumerate(input_string):
                        dna_sequences.append(string)
                else:
                    raise ValueError("There is no digital data input here!")

                original_dna_sequences = copy.deepcopy(dna_sequences)

                results = self.coding_scheme.carbon_to_silicon(dna_sequences)
                self.records["decoding runtime"] = round(results["t"], 3)

                bit_segments = results["bit"]
                bit_size = results["s"]

                if not bit_segments:
                    self.records["error rate"] = "100.00%"
                    return {"bit": None, "dna": original_dna_sequences}

                if self.error_correction is not None:
                    verified_data = self.error_correction.remove(bit_segments)
                    bit_segments = verified_data["bit"]
                    self.records["error rate"] = str(round(verified_data["e_r"] * 100, 2)) + "%"
                    self.records["error indices"] = str(verified_data["e_i"]).replace(", ", "-") \
                        if verified_data["e_i"] != [] else None
                    self.records["error bit segments"] = str(verified_data["e_bit"]).replace(", ", "-") \
                        if verified_data["e_bit"] != [] else None
                else:
                    self.records["error rate"] = None
                    self.records["error indices"] = None
                    self.records["error bit segments"] = None

                if not bit_segments:
                    return {"bit": None, "dna": original_dna_sequences}

                if "index" in info and info["index"]:
                    if index_length is not None:
                        indices, bit_segments = indexer.divide_all(bit_segments, index_length, self.need_logs)
                    else:
                        indices, bit_segments = indexer.divide_all(bit_segments, None, self.need_logs)

                    bit_segments = indexer.sort_order(indices, bit_segments, self.need_logs)

                if "output_path" in info:
                    data_handle.write_bits_to_file(info["output_path"], bit_segments, bit_size, self.need_logs)
                elif "output_string" in info:
                    string = data_handle.write_bits_to_str(bit_segments, bit_size, self.need_logs)
                    if self.need_logs:
                        print(string)

                return {"bit": bit_segments, "size": bit_size, "dna": original_dna_sequences}
            else:
                raise ValueError("Unknown parameter \"direction\", please use \"t_c\" or \"t_s\".")
        else:
            raise ValueError("Unknown parameter \"direction\", please use \"t_c\" or \"t_s\".")
