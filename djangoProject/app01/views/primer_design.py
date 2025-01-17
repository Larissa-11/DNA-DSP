import json
import re
from threading import Thread

from Bio.Seq import Seq
from django.http import JsonResponse
from django.views import View

from app01.Primer_design.collisions import primer_payload_collision
from app01.Primer_design.primer_design import design_primers
from app01.Primer_design.remove_primers import remove_primers_from_payloads
from app01.models import TbResult
from app01.views.algorithm_evaluation import thread_pool
from app01.views.algorithm_evaluation import notification_queue
from app01.views.errorInfoList import ERROR_INFO_LIST
from djangoProject.tools.webSocket import send_notification


class primer_design(View):

    def post(self, request):
        try:

            params = request.POST or json.loads(request.body.decode())
            keys = params.get("keys")
            user_id = request.user.id
            result = {}
            msg = "OK"
            if keys == "1":
                # 处理参数类型转换
                p = {
                    'length': int(params.get('primer_length')),
                    'gc_content_min': float(params.get('gc_min')),
                    'gc_content_max': float(params.get('gc_max')),
                    'homo_max_len': int(params.get('homopolymer_length')),
                    'Tm_min': float(params.get('mt_min')),
                    'Tm_max': float(params.get('mt_max')),
                    'hamming_distance': int(params.get('homopolymer_distance')),
                    'inter_complementarity': int(params.get('inter_complementarity')),
                    'intra_complementarity': int(params.get('Intra_complementarity')),
                    'number': int(params.get('primer_number'))
                }
                print(p)

                data = {
                    "p": p,

                }
                if user_id:
                    data['user_id'] = user_id
                    thread_pool.submit(primer_design_func, **data)
                else:
                    result, msg = primer_design_func_01(**data)

            elif keys == "3":
                kwy_file = request.FILES['key_file']
                p = {
                    'forward_primer_length': int(params.get('forward_primer1', 20)),
                    'reverse_primer_length': int(params.get('reverse_primer2', 20)),
                    'payload_data': kwy_file.read().decode('utf-8'),
                }
                data = {
                    "p": p,
                    "content_type": kwy_file.content_type,
                    "file_name": kwy_file.name,
                }
                if user_id:
                    data['user_id'] = user_id
                    thread_pool.submit(remove_primers_func, **data)
                else:
                    result, msg = remove_primers_func_01(**data)
            elif keys == "2":
                kwy_file = request.FILES['key_file']
                data = {
                    'forward_primer': params.get('forward_primer'),
                    'reverse_primer': params.get('reverse_primer'),
                    'payloads': kwy_file.read().decode('utf-8'),

                    "content_type": kwy_file.content_type,
                    "file_name": kwy_file.name,
                }
                if user_id:
                    data['user_id'] = user_id
                    thread_pool.submit(collision_func, **data)
                else:
                    result, msg = collision_func_01(**data)
            return JsonResponse({"code": 200, "data": json.dumps(result), "msg": msg})
        except Exception as e:
            print("ERROR:", e)
            return JsonResponse({"code": 200, "data": None, "msg": "Error"})


def primer_design_func(p, user_id):
    try:
        primers = design_primers(p)
        results = ""
        for i, primer in enumerate(primers):
            results += primer + "\n"

        data = [
            {
                "fileData": results,
                "fileName": "primer_design.txt",
                "fileType": "txt",
            }
        ]

        tmp_model = TbResult.objects.create(result=json.dumps(data),
                                            user_id=user_id,
                                            dataType="Primer_design")

        send_notification(tmp_model.id, "Primer_design", "Primer_design", data=json.dumps(data))
    except Exception as e:
        print("error:", e)
        error_info = str(e)
        key_name = "PrimerDesignError"
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


def primer_design_func_01(p):
    primers = design_primers(p)
    results = ""
    for i, primer in enumerate(primers):
        results += primer + "\n"

    data = [
        {
            "fileData": results,
            "fileName": "primer_design.txt",
            "fileType": "txt",
        }
    ]
    return data, "OK"


def remove_primers_func(user_id, content_type, file_name, p):
    try:
        results = remove_primers_from_payloads(**p)

        data = [
            {
                "fileData": results,
                "fileName": file_name.split(".")[0] + "_Remove Primer" + ".txt",
                "fileType": "txt",
            }
        ]

        tmp_model = TbResult.objects.create(result=json.dumps(data),
                                            user_id=user_id,
                                            dataType="Remove primers")

        send_notification(tmp_model.id, "Remove primers", "Remove primers", data=json.dumps(data))
    except Exception as e:
        print("error:", e)
        error_info = str(e)
        key_name = "RemovePrimersError"
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


def remove_primers_func_01(content_type, file_name, p):
    try:
        results = remove_primers_from_payloads(**p)

        data = [
            {
                "fileData": results,
                "fileName": file_name.split(".")[0] + "_Remove Primer" + ".txt",
                "fileType": "txt",
            }
        ]
        return data, None
    except Exception as e:
        print("error:", e)


def collision_func(user_id, content_type, file_name, payloads, forward_primer, reverse_primer):
    try:
        # 判断一下forward_primer和reverse_primer的长度，如果长度都小于12，那直接返回一个消息，提示，无需进行碰撞检测。
        # Convert the reverse primer to its reverse complement
        # reverse_primer_rc = str(Seq(reverse_primer).reverse_complement())
        payloads = re.sub(r'\n+$', '', payloads)
        payloads = payloads.split("\n")
        # Check for collisions
        print("Checking forward primer for collisions...")
        forward_collision = primer_payload_collision(forward_primer, payloads)

        print("Checking reverse primer for collisions...")
        reverse_collision = primer_payload_collision(reverse_primer, payloads)
        # forward_collision = None
        # reverse_collision = 1
        if forward_collision or reverse_collision:
            tmp_model = TbResult.objects.create(result="",
                                                user_id=user_id,
                                                dataType="Collision")

            send_notification(tmp_model.id, "Collision", "Collision", data="")

        else:
            modified_payloads = ""
            for payload in payloads:
                modified_payloads += forward_primer + payload + reverse_primer + "\n"
            # modified_payloads = [forward_primer + payload + reverse_primer for payload in payloads]
            # modified_payloads
            data = [
                {
                    "fileData": modified_payloads,
                    "fileName": file_name.split(".")[0] + "_Primer" + ".txt",
                    "fileType": "txt",
                }
            ]
            tmp_model = TbResult.objects.create(result=json.dumps(data),
                                                user_id=user_id,
                                                dataType="Add primers")

            send_notification(tmp_model.id, "Add primers", "Add primers", data=json.dumps(data))
    except Exception as e:
        print("error:", e)
        error_info = str(e)
        key_name = "AddPrimersError"
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


def collision_func_01(content_type, file_name, payloads, forward_primer, reverse_primer):
    try:
        # 判断一下forward_primer和reverse_primer的长度，如果长度都小于12，那直接返回一个消息，提示，无需进行碰撞检测。
        # Convert the reverse primer to its reverse complement
        reverse_primer_rc = str(Seq(reverse_primer).reverse_complement())
        payloads = re.sub(r'\n+$', '', payloads)
        payloads = payloads.split("\n")
        # Check for collisions
        print("Checking forward primer for collisions...")
        forward_collision = primer_payload_collision(forward_primer, payloads)

        print("Checking reverse primer for collisions...")
        reverse_collision = primer_payload_collision(reverse_primer_rc, payloads)
        # forward_collision = None
        # reverse_collision = 1
        if forward_collision or reverse_collision:
            return None, "Collision"
        else:
            modified_payloads = ""
            for payload in payloads:
                modified_payloads += forward_primer + payload + reverse_primer + "\n"
            # modified_payloads = [forward_primer + payload + reverse_primer for payload in payloads]
            # modified_payloads
            data = [
                {
                    "fileData": modified_payloads,
                    "fileName": file_name.split(".")[0] + "_Primer" + ".txt",
                    "fileType": "txt",
                }
            ]
            return data, "OK"
    except Exception as e:
        print("error:", e)
        return None, f"{e}"