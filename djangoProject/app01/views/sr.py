import json
from threading import Thread

from django.utils.functional import classproperty
from django.views import View
from django.http import JsonResponse
import time
import os

from app01.models import TbResult
from app01.sequence_reconstruction.clover.clover_clustering import clover_method
from app01.sequence_reconstruction.fuzzy_clutering.HSFC_method import HSFC_method
from app01.sequence_reconstruction.lsh.lsh import LSH_method
from app01.sequence_reconstruction.msa.msa import process_byte_stream, clusters_msa, generate_candidates, \
    format_sequences_as_fasta
from app01.views.algorithm_evaluation import thread_pool
from app01.views.algorithm_evaluation import notification_queue
from app01.views.coding_scheme import decode_dask
from app01.views.errorInfoList import ERROR_INFO_LIST
from djangoProject.settings import BASE_DIR
from djangoProject.tools.webSocket import send_notification


class ClustringView(View):
    @classproperty
    def clustering_method_map(cls) -> dict:
        return {
            "clover_method": cls.clover_process,
            "lsh_method": cls.lsh_process,
            "HSFC_method": cls.HSFC_process,
        }

    @classproperty
    def MSA_method_map(cls) -> dict:
        return {
            "MSA_method": cls.muscle_process,
        }

    @staticmethod
    def clover_process(request):

        data = request.POST or json.loads(request.body.decode())
        read_len = data.get("read_len")
        end_tree_len = data.get("end_tree_len")
        other_tree_len = data.get("other_tree_len")
        Vertical_drift = data.get("Vertical_drift")
        Horizontal_drift = data.get("Horizontal_drift")
        Cluster_size_threshold = data.get("Cluster_size_threshold")
        h_index_nums = data.get("h_index_nums")
        e_index_nums = data.get("e_index_nums")
        read_len_min = data.get("read_len_min")
        try:
            read_len = int(read_len)
            end_tree_len = int(end_tree_len)
            other_tree_len = int(other_tree_len)
            Vertical_drift = int(Vertical_drift)
            Horizontal_drift = int(Horizontal_drift)
            Cluster_size_threshold = int(Cluster_size_threshold)
            h_index_nums = int(h_index_nums)
            e_index_nums = int(e_index_nums)
            read_len_min = int(read_len_min)
        except:
            return JsonResponse({
                "code": 400,
                "data": None,
                "msg": "Missing read_len parameter"
            })
        file = request.FILES.get('file')
        byte_stream = file.read()
        input_path = f"input{time.time()}.txt"
        input_path = os.path.join(BASE_DIR, 'templates', input_path)

        key = {
            "input_path": input_path,
            "byte_stream": byte_stream,
            "read_len": read_len,
            "end_tree_len": end_tree_len,
            "other_tree_len": other_tree_len,
            "Vertical_drift": Vertical_drift,
            "Horizontal_drift": Horizontal_drift,
            "Cluster_size_threshold": Cluster_size_threshold,
            "h_index_nums": h_index_nums,
            "e_index_nums": e_index_nums,
            "read_len_min": read_len_min,
            "user_id": request.user.id,
            "file_name": file.name,
        }

        def task(input_path, byte_stream, read_len, end_tree_len, other_tree_len, Vertical_drift,
                 Horizontal_drift, Cluster_size_threshold, h_index_nums, e_index_nums, read_len_min, user_id,
                 file_name, ):
            try:
                with open(input_path, 'wb') as f:
                    f.write(byte_stream)
                results_string_noseq, results_list = clover_method(input_path, "output",
                                                                   read_len, end_tree_len, other_tree_len,
                                                                   Vertical_drift,
                                                                   Horizontal_drift, Cluster_size_threshold,
                                                                   h_index_nums,
                                                                   e_index_nums, read_len_min)

                data = [
                    {
                        "fileData": results_list,
                        "fileName": file_name.split(".")[0] + "_Clover",
                        "fileType": file_name.split(".")[-1],
                    },
                    {
                        "fileData": results_string_noseq,
                        "fileName": file_name.split(".")[0] + "_Clover_Noseq",
                        "fileType": file_name.split(".")[-1],
                    }
                ]
                tmp_model = TbResult.objects.create(result=json.dumps(data),
                                                    user_id=user_id,
                                                    dataType="Clustering_clover")
                send_notification(tmp_model.id, "Clustering_clover", "Clustering_clover", data=json.dumps(data))
            except Exception as e:
                error_info = str(e)
                print("ClusteringCloverError", error_info)
                key_name = "ClusteringCloverError"
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

        # thread = Thread(target=task, kwargs=key)
        # # 启动线程
        # thread.start()
        thread_pool.submit(task, **key)
        return JsonResponse({
            "code": 200,
            "data": None,
            "msg": "Success",
            "content_type": "",
            "suffix_name": ""
        })

    @staticmethod
    def lsh_process(request):
        file = request.FILES.get('file')
        byte_stream = file.read()
        key = {
            "file_name": file.name,
            "user_id": request.user.id,
            "byte_stream": byte_stream,
        }

        def task(byte_stream, file_name, user_id):
            try:
                sequences = [_.strip() for _ in byte_stream.decode().split("\n") if _]
                results_string_noseq, results_string = LSH_method(sequences)

                data = [
                    {
                        "fileData": results_string,
                        "fileName": file_name.split(".")[0] + "_LSH",
                        "fileType": file_name.split(".")[-1],
                    },
                    {
                        "fileData": results_string_noseq,
                        "fileName": file_name.split(".")[0] + "_LSH_Noseq",
                        "fileType": file_name.split(".")[-1],
                    }
                ]
                tmp_model = TbResult.objects.create(result=json.dumps(data),
                                                    user_id=user_id,
                                                    dataType="Clustering_lsh")
                send_notification(tmp_model.id, "Clustering_lsh", "Clustering_lsh", data=json.dumps(data))
            except Exception as e:
                error_info = str(e)
                key_name = "ClusteringLshError"
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

        # thread = Thread(target=task, kwargs=key)
        # # 启动线程
        # thread.start()
        thread_pool.submit(task, **key)
        # results_string_noseq  = '\n'.join(results_string_noseq)
        # results_string = '\n'.join(results_string)
        return JsonResponse({
            "code": 200,
            "data": None,
            "msg": "Success",
            "content_type": file.content_type,
            "suffix_name": file.name.split(".")[-1]
        })

    @staticmethod
    def HSFC_process(request):
        data = request.POST or json.loads(request.body.decode())
        correct_length = data.get("correct_length")
        try:
            correct_length = int(correct_length)
        except:
            return JsonResponse({
                "code": 400,
                "data": None,
                "msg": "Missing read_len parameter"
            })
        file = request.FILES.get('file')
        byte_stream = file.read()
        key = {
            "file_name": file.name,
            "user_id": request.user.id,
            "byte_stream": byte_stream,
            "correct_length": correct_length,
        }

        def task(byte_stream, file_name, user_id, correct_length):
            try:
                sequences = [_.strip() for _ in byte_stream.decode().split("\r\n") if _]
                results_string_noseq, results_string = HSFC_method(sequences, correct_length)
                data = [
                    {
                        "fileData": results_string,
                        "fileName": file_name.split(".")[0] + "_HSFC",
                        "fileType": file_name.split(".")[-1],
                    },
                    {
                        "fileData": results_string_noseq,
                        "fileName": file_name.split(".")[0] + "_HSFC_Noseq",
                        "fileType": file_name.split(".")[-1],
                    }
                ]
                tmp_model = TbResult.objects.create(result=json.dumps(data),
                                                    user_id=user_id,
                                                    dataType="Clustering_HSFC")
                send_notification(tmp_model.id, "Clustering_HSFC", "Clustering_HSFC", data=json.dumps(data))
            except Exception as e:
                error_info = str(e)
                key_name = "ClusteringHSFCError"
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

        # thread = Thread(target=task, kwargs=key)
        # # 启动线程
        # thread.start()
        thread_pool.submit(task, **key)
        # results_string_noseq = '\n'.join(results_string_noseq)
        # results_string = '\n'.join(results_string)
        return JsonResponse({
            "code": 200,
            "data": None,
            "msg": "Success",
            "content_type": file.content_type,
            "suffix_name": file.name.split(".")[-1]
        })

    @staticmethod
    def muscle_process(request):
        file = request.FILES.get('file')
        byte_stream = file.read()
        key = {
            "file_name": file.name,
            "user_id": request.user.id,
            "byte_stream": byte_stream,
        }

        def task(byte_stream, file_name, user_id, ):
            try:
                reads, cluster_results = process_byte_stream(byte_stream)
                msa_results = clusters_msa(cluster_results, reads)
                candidates_results = generate_candidates(msa_results)
                results_string = format_sequences_as_fasta(candidates_results)

                data = [
                    {
                        "fileData": results_string,
                        "fileName": file_name.split(".")[0] + "_MSA",
                        "fileType": file_name.split(".")[-1],
                    }
                ]
                tmp_model = TbResult.objects.create(result=json.dumps(data),
                                                    user_id=user_id,
                                                    dataType="MSA_Muscle")
                send_notification(tmp_model.id, "MSA_Muscle", "MSA_Muscle", data=json.dumps(data))
            except Exception as e:
                print(str(e))
                error_info = str(e)
                key_name = "MSA_Muscle_Error"
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

        # thread = Thread(target=task, kwargs=key)
        # # 启动线程
        # thread.start()
        thread_pool.submit(task, **key)
        return JsonResponse({
            "code": 200,
            "data": None,
            "msg": "Success",
            "content_type": file.content_type,
            "suffix_name": file.name.split(".")[-1]
        })

    def post(self, request):
        data = request.POST or json.loads(request.body.decode())
        # 这个应该可以直接去掉，在前端不是判断了吗？
        if not request.FILES.get('file'):
            return JsonResponse({
                "code": 400,
                "data": None,
                "msg": "Inexistence"
            })
        clustering_value = request.POST.get('Clustering')
        MSA_value = request.POST.get('MSA')

        decode_flag = request.POST.get('decode_flag')
        if clustering_value == 'true':
            clustering_method = data.get("clustering_method")
            process_method = self.clustering_method_map.get(clustering_method)
        elif MSA_value == 'true':
            msa_method = data.get("MSA_method")
            process_method = self.MSA_method_map.get(msa_method)
        elif decode_flag == 'true':
            process_method = ""
            user_id = request.user.id
            # user_id = 1
            data = request.POST or json.loads(request.body.decode())
            file = request.FILES.get('file')

            data = {
                "errorCorrection": data.get('errorCorrection'),
                "codingScheme": data.get('codingScheme'),
                "neededIndice": data.get('neededIndice'),
                "segmentLength": data.get('segmentLength', 0),
            }

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
            return JsonResponse({
                "code": 200,
                "data": None,
                "msg": "OK"
            })
        else:
            return JsonResponse({
                "code": 400,
                "data": None,
                "msg": "Please select a method"
            })
        if not process_method:
            return JsonResponse({
                "code": 400,
                "data": None,
                "msg": "Please select a method"
            })
        return process_method(request)
