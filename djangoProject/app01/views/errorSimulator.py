import json
from threading import Thread

from django.http import HttpResponse, JsonResponse
from django.views import View

from app01.ErrorSimulator.error_simulator import parse_rates_dictionary, Simulator
from app01.models import TbResult
from app01.views.algorithm_evaluation import thread_pool
from app01.views.algorithm_evaluation import notification_queue
from app01.views.errorInfoList import ERROR_INFO_LIST
from djangoProject.tools.webSocket import send_notification


class ErrorSimulator(View):
    def post(self, request):
        # try:
        user = request.user
        user_id = user.id
        # user_id = 1
        total_error_rates_example = {
            'd': float(request.POST.get("total_error_rates_example[One_Base_Deletion]")),
            'ld': float(request.POST.get("total_error_rates_example[Long_Deletion]")),
            'i': float(request.POST.get("total_error_rates_example[Insertion]")) ,
            's': float(request.POST.get("total_error_rates_example[Substitution]")) }
        base_error_rates_example = {'A':
                                        {'s': float(request.POST.get('key_A[substitution]')) ,
                                         'i': float(request.POST.get("key_A[insertion]")) ,
                                         'pi': float(request.POST.get("key_A[pre_insertion]")),
                                         'd': float(request.POST.get("key_A[base_deletion]")) ,
                                         'ld': float(request.POST.get("key_A[long_deletion]"))},
                                    'C':
                                        {'s': float(request.POST.get("key_C[substitution]")),
                                         'i': float(request.POST.get("key_C[insertion]")),
                                         'pi': float(request.POST.get("key_C[pre_insertion]")),
                                         'd': float(request.POST.get("key_C[base_deletion]")),
                                         'ld': float(request.POST.get("key_C[long_deletion]"))},
                                    'T':
                                        {'s': float(request.POST.get("key_T[substitution]")),
                                         'i': float(request.POST.get("key_T[insertion]")),
                                         'pi': float(request.POST.get("key_T[pre_insertion]")),
                                         'd': float(request.POST.get("key_T[base_deletion]")) ,
                                         'ld': float(request.POST.get("key_T[long_deletion]"))},
                                    'G':
                                        {'s': float(request.POST.get("key_G[substitution]")),
                                         'i': float(request.POST.get("key_G[insertion]")) ,
                                         'pi': float(request.POST.get("key_G[pre_insertion]")),
                                         'd': float(request.POST.get("key_G[base_deletion]")),
                                         'ld': float(request.POST.get("key_G[long_deletion]"))},
                                    }
        parse_rates_dictionary(total_error_rates_example)
        parse_rates_dictionary(base_error_rates_example)
        min_copies = float(request.POST.get("Min"))
        max_copies = float(request.POST.get("Max"))
        #
        input_path_example = request.FILES.getlist('input_path_example[0]')
        #
        #
        if input_path_example is None:
            return

        input_path_example = input_path_example[0]

        key = {
            "base_error_rates_example": base_error_rates_example,
            "input_path_example": input_path_example,
            "max_copies": max_copies,
            "min_copies": min_copies,
            "total_error_rates_example": total_error_rates_example,

            "file_name": input_path_example.name.split('.')[0],
            "file_type":  input_path_example.name.split('.')[1]
        }
        # self.asyc_task(base_error_rates_example, input_path_example, max_copies, min_copies, total_error_rates_example)
        # thread = Thread(target=self.asyc_task, kwargs=key)
        # # 启动线程
        # thread.start()
        if user_id:
            key['user_id'] = user_id
            thread_pool.submit(self.asyc_task, **key)
            return JsonResponse({
                "code": 200,
                "data": "",
                "msg": "",

            })
        else:
            result, msg = self.asyc_task_01(**key)
            return JsonResponse({
                "code": 200,
                "data": json.dumps(result),
                "msg": msg,

            })

    def asyc_task(self, base_error_rates_example, input_path_example, max_copies, min_copies,
                  total_error_rates_example, user_id, file_name, file_type):
        try:
            sim = Simulator(total_error_rates_example, base_error_rates_example, min_copies, max_copies,
                            input_path_example)
            evyat_str, shuffled_str = sim.simulate_errors()


            tmp_list = shuffled_str.split("\n")

            file_one = ""
            file_two = ""

            j = 0
            k = 1
            for i in tmp_list:
                file_one += ">Seq" + str(j) + "\n" + i + "\n"

                file_two += str(k) + " " + i + "\n"

                j+=1
                k+=1


            data = [{
                    "fileData": file_one,
                    "fileName": file_name + "_ErrorSimulator." + "txt",
                    "fileType": "txt",
                },{
                    "fileData": file_two,
                    "fileName": file_name + "_ErrorSimulator_Noseq." + "txt",
                    "fileType": "txt",
                },]

            tmp_model = TbResult.objects.create(result=json.dumps(data),
                                                user_id=user_id,
                                                dataType="ErrorSimulator")

            send_notification(tmp_model.id, "ErrorSimulator", "ErrorSimulator", data=json.dumps(data))
        except Exception as e:
            print("error:", e)
            error_info = str(e)
            if error_info not in ERROR_INFO_LIST:
                error_info = "ErrorSimulatorError"
            tmp_model = TbResult.objects.create(result="",
                                                user_id=user_id,
                                                error=error_info,
                                                dataType="ErrorSimulatorError")

            # 将消息放入队列
            notification_queue.put(
                {"id": tmp_model.id, "error": error_info, "message": "ErrorSimulatorError", "data": None,
                 "type_tmp": "ErrorSimulatorError"})

    def asyc_task_01(self, base_error_rates_example, input_path_example, max_copies, min_copies,
                  total_error_rates_example, file_name, file_type):
        try:
            sim = Simulator(total_error_rates_example, base_error_rates_example, min_copies, max_copies,
                            input_path_example)
            evyat_str, shuffled_str = sim.simulate_errors()


            tmp_list = shuffled_str.split("\n")

            file_one = ""
            file_two = ""

            j = 0
            k = 1
            for i in tmp_list:
                file_one += ">Seq" + str(j) + "\n" + i + "\n"

                file_two += str(k) + " " + i + "\n"

                j+=1
                k+=1


            data = [{
                    "fileData": file_one,
                    "fileName": file_name + "_ErrorSimulator." + "txt",
                    "fileType": "txt",
                },{
                    "fileData": file_two,
                    "fileName": file_name + "_ErrorSimulator_Noseq." + "txt",
                    "fileType": "txt",
                },]
            return data, "OK"
        except Exception as e:
            print("error:", e)
            return None, f"{e}"

# except Exception as e:
#     print(e)
#     return HttpResponse({
#         "code": 999999,
#         "msg": "error",
#     })
