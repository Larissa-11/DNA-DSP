import base64
import json
import queue
import struct
import threading
import time
import uuid
from threading import Thread

import numpy as np
from django.db import IntegrityError

from Chamaeleo.methods.default import BaseCodingAlgorithm
from Chamaeleo.methods.ecc import Hamming, ReedSolomon
from Chamaeleo.methods.fixed import Church, Goldman, Grass, Blawat
from Chamaeleo.methods.flowed import DNAFountain, YinYangCode
from app01.models import DnaInfo, TbResult
from app01.views.algorithm_evaluation import notification_queue
from app01.views.errorInfoList import ERROR_INFO_LIST
from djangoProject.tools.webSocket import send_notification

from utils.auth import LoginRequiredJsonMixin
from utils.my_chamaeleo.pipelines import MyTranscodePipeline
from django.http import JsonResponse
from django.views import View
from django.core.handlers.wsgi import WSGIRequest


def success_res(data=None, msg: str = "success", code=200, **kwargs):
    return JsonResponse({
        "code": code,
        "data": data,
        "msg": msg,
        **kwargs
    })


def error_res(msg: str = "error", code=400, **kwargs):
    return JsonResponse({
        "code": code,
        "data": None,
        "msg": msg,
        **kwargs
    })


def handle_combined_transcode(direction: str, data, byte_stream, user, file_name, content_type) -> dict:
    try:
        print("handle_combined_transcode Start")
        error_correction_str = data.get('errorCorrection')
        coding_scheme_str = data.get('codingScheme')
        needed_indice = data.get('neededIndice')
        segment_length = data.get('segmentLength', 0)
        # index_length = data.get('indexLength')
        index_length = None
        Bits = data.get('Bits', 3)
        if type(Bits) is str:
            Bits = int(Bits)
        # if index_length.lower() == "nan" else index_length

        if None in [error_correction_str, coding_scheme_str, needed_indice, segment_length]:
            raise ValueError("Incomplete parameters")

        try:
            if not isinstance(needed_indice, bool):
                needed_indice = needed_indice.lower() == "true"
            segment_length = int(segment_length)
            index_length = int(index_length) if index_length else None
            index_length = None if index_length == 0 else index_length
        except:
            raise ValueError("Parameter type error")

        # 各个编码方案映射
        CODING_SCHEME_MAP = {
            "Base": BaseCodingAlgorithm(),
            "Church": Church(),
            "Goldman": Goldman(),
            "Grass": Grass(),
            "Blawat": Blawat(),
            "DNAFountain": DNAFountain(redundancy=0.5),
            "YinYangCode": YinYangCode(),
        }

        ERROR_CORRECTION_MAP = {
            "None": None,
            "Hamming": Hamming(),
            "ReedSolomon": ReedSolomon(Bits)
            # "ReedSolomon": ReedSolomon(6)
        }

        data_dict = {}
        coding_scheme = CODING_SCHEME_MAP.get(coding_scheme_str)
        prefix_name = file_name.split(".")[0]
        suffix_name = file_name.split(".")[-1]
        # 解码需要赋值元数据
        if direction == "t_s":
            try:
                info = DnaInfo.objects.get(file_name=prefix_name, user=user)
            except DnaInfo.DoesNotExist:
                raise ValueError("This dna file is not encoded by this system！")
            data_dict["content_type"] = info.content_type
            data_dict["suffix_name"] = info.suffix_name
            coding_scheme.bit_size = info.bit_size
            coding_scheme.segment_length = info.segment_length
            coding_scheme.index_length = info.index_length
            segment_length = info.prime_segment_length
            total_count = info.total_count
            ERROR_CORRECTION_MAP["ReedSolomon"] = ReedSolomon(info.bits)
            if coding_scheme_str == "DNAFountain":
                coding_scheme.decode_packets = total_count
            else:
                coding_scheme.total_count = total_count

        error_correction = ERROR_CORRECTION_MAP.get(error_correction_str)
        pipeline = MyTranscodePipeline(
            coding_scheme=coding_scheme,
            error_correction=error_correction
        )
        try:
            data_dict.update(
                pipeline.transcode(
                    direction=direction,
                    input_byte=byte_stream,
                    segment_length=segment_length,
                    index_length=index_length,
                    index=needed_indice
                )
            )
        except AttributeError as e:
            raise e

        # 如果是编码把一些原数据保存到数据库
        if direction == "t_c":
            # 文件类型
            file_name = f'{prefix_name}_{coding_scheme_str}{(int(time.time()))}'
            try:
                DnaInfo.objects.create(
                    user=user,
                    content_type=content_type,
                    file_name=file_name,
                    suffix_name=suffix_name,
                    bits=Bits,
                    prime_segment_length=segment_length,
                    **data_dict["info"]
                )
            except IntegrityError:
                raise ValueError("The file already exists, please recode it!")
            data_dict["file_name"] = file_name
        print("handle_combined_transcode End")
        return data_dict
    except Exception as e:
        print("errpr", e)
        raise e


class EncodeView(View):
    def post(self, request):
        try:

            user_id = request.user.id
            # user_id = 1
            data = request.POST or json.loads(request.body.decode())
            file = request.FILES.get('file')
            key = {
                "user_id": user_id,
                "user": request.user,
                "data": data,
                "byte_stream": file.read(),
                "content_type": file.content_type,
                "file_name": file.name,
            }

            thread = Thread(target=task, kwargs=key)
            # 启动线程
            thread.start()
            return success_res()
        except Exception as e:  # noqa
            return error_res(str(e))


def task(user_id, data, byte_stream, user, content_type, file_name):
    error_info = ""
    try:
        print("EncodeView Start")
        scheme_data_dict = handle_combined_transcode(direction="t_c", data=data,
                                                     byte_stream=byte_stream, user=user,
                                                     content_type=content_type, file_name=file_name
                                                     )
        encoded_data = scheme_data_dict['dna']
        result = "\n".join("".join(i) for i in encoded_data)
        data1 = [
            {
                "fileData": result,
                "fileName": scheme_data_dict["file_name"],
                "fileType": "txt",
            }
        ]
        tmp_model = TbResult.objects.create(result=json.dumps(data1),
                                            user_id=user_id,
                                            dataType="Encode")
        send_notification(tmp_model.id, "Encode", "Encode", data=json.dumps(data1))
        print("EncodeView End")
    except Exception as e:
        print(e)
        key_name = "EncodeError"
        error_info = str(e)
        if error_info not in ERROR_INFO_LIST:
            error_info = "EncodeError"
        tmp_model = TbResult.objects.create(result="",
                                            user_id=user_id,
                                            error=str(error_info),
                                            dataType=key_name)

        # 将消息放入队列
        notification_queue.put(
            {"id": tmp_model.id, "error": error_info, "message": key_name, "data": None,
             "type_tmp": key_name})
    # finally:
    #     tmp_model = TbResult.objects.create(result="",
    #                                         user_id=user_id,
    #                                         error=str(error_info),
    #                                         dataType="Encode_Error")
    #     send_notification(tmp_model.id, "Encode", "Encode", error=str(error_info), data="")


class DecodeView(View):
    def post(self, request):
        try:
            user_id = request.user.id
            # user_id = 1
            data = request.POST or json.loads(request.body.decode())
            file = request.FILES.get('file')

            key = {
                "user_id": user_id,
                "user": request.user,
                "data": data,
                "byte_stream": file.read(),
                "content_type": file.content_type,
                "file_name": file.name.split("_")[0],
                "time_Str": file.name.split("_")[1].split(".")[0],
                "file_type": file.name.split(".")[-1],
            }

            thread = Thread(target=decode_dask, kwargs=key)
            # 启动线程
            thread.start()
            return success_res()

        except Exception as e:
            return error_res(str(e))


def decode_dask(user_id, data, byte_stream, user, content_type, file_name, time_Str, file_type):
    try:
        file_name_tmp = file_name + "_" + time_Str + "." + file_type
        scheme_data_dict = handle_combined_transcode(direction="t_s", data=data,
                                                     byte_stream=byte_stream, user=user,
                                                     content_type=content_type, file_name=file_name_tmp
                                                     )
        bit = scheme_data_dict["bit"]
        size = scheme_data_dict.get("size")
        if bit is None:
            return error_res("decode error")

        # 将二进制矩阵转换为字节串
        values = []
        byte_string = b''
        matrix = np.array(bit).reshape(-1)
        for position in range(0, size, 8):
            byte_string += struct.pack("B", int("".join(list(map(str, matrix[position: position + 8]))), 2))
        # for index in range(0, len(matrix), 8):
        #     values.append(int("".join(list(map(str, matrix[index: index + 8]))), 2))
        # base64_str=str(bytes(values), encoding="utf8")
        # 转为base64编码
        base64_str = base64.b64encode(byte_string).decode('utf-8')
        content_type = scheme_data_dict["content_type"]
        suffix_name = scheme_data_dict["suffix_name"]
        data = [
            {
                "fileData": base64_str,
                "fileName": file_name,
                "fileType": content_type,
            }
        ]
        tmp_model = TbResult.objects.create(result=json.dumps(data),
                                            user_id=user_id,
                                            dataType="Decode")
        send_notification(tmp_model.id, "Decode", "Decode", data=json.dumps(data))
    except Exception as e:
        print("error:", e)
        error_info = str(e)
        key_name = "Decode_Error"
        if error_info not in ERROR_INFO_LIST:
            error_info = key_name
        tmp_model = TbResult.objects.create(result="",
                                            user_id=user_id,
                                            error=error_info,
                                            dataType=key_name)

        # 将消息放入队列
        notification_queue.put(
            {"id": tmp_model.id, "error": error_info, "message": key_name, "data": None,
             "type_tmp": key_name})
