import json
import queue
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from threading import Thread, current_thread

from django.http import JsonResponse
from django.views import View

from app01.encoding_analysis.Net_hom_motifs import CalculateNetHomopolymerMotifs
from app01.encoding_analysis.GC_content_distribution import GcContentDistribution
from app01.encoding_analysis.MFE import Mfe
from djangoProject.tools.webSocket import send_notification
from ..models import TbResult

thread_pool = ThreadPoolExecutor(max_workers=20)
# 创建通知队列
notification_queue = queue.Queue()
def notification_sender():
    while True:
        try:
            msg = notification_queue.get(timeout=1)  # 等待新消息
            send_notification(msg['id'], msg['message'], msg['type_tmp'], error=msg['error'], data=msg['data'])
        except queue.Empty:
            continue
# 启动通知发送线程
threading.Thread(target=notification_sender, daemon=True).start()


class AlgorithmEvaluationView(View):
    @staticmethod
    def get_dna_sequence(file):
        byte_stream = file.read()
        lines = byte_stream.decode().split('\n')
        dna_sequences = [list(line.strip()) for line in lines if line.strip()]
        return dna_sequences

    def post(self, request):
        original_file = request.FILES.getlist('original_file')
        # original_file = [file.name for file in original_file]
        file_list = request.FILES.getlist('file')
        file_name_list = [file.name for file in file_list]
        bit_size = CalculateNetHomopolymerMotifs.read_bits_from_file(original_file[0])
        dna_sequences_list = [self.get_dna_sequence(file) for file in file_list]

        net_info_density = CalculateNetHomopolymerMotifs.calculate_net_information_density(bit_size, dna_sequences_list)
        max_lengths = CalculateNetHomopolymerMotifs.find_longest_homopolymer(dna_sequences_list)
        motifs_count = CalculateNetHomopolymerMotifs.motifs(dna_sequences_list)
        gc_content_distribution = GcContentDistribution()
        gc_content_result = gc_content_distribution.gc_content_distribution(dna_sequences_list)
        local_gc_content_result = gc_content_distribution.local_gc_content_distribution(dna_sequences_list)
        mfe = Mfe()
        mfe_result = mfe.mfe_calculation(dna_sequences_list)
        data = {"xAxis": file_name_list,
                "yAxis": {"net_info_density": net_info_density,
                          "max_lengths": max_lengths,
                          "motifs_count": motifs_count,
                          "gc_content_result": gc_content_result,
                          "mfe_result": mfe_result,
                          "local_gc_content_result": local_gc_content_result
                          }}
        return JsonResponse({
            "code": 200,
            "data": data,
            "msg": "msg",
        })


from ..tasks import test03, test02


class Test(View):
    @staticmethod
    def get_dna_sequence(file):
        byte_stream = file.read()
        lines = byte_stream.decode().split('\n')
        dna_sequences = [list(line.strip()) for line in lines if line.strip()]
        return dna_sequences

    def post(self, request):
        original_file = request.FILES.getlist('original_file')
        # original_file = [file.name for file in original_file]
        file_list = request.FILES.getlist('file')
        file_name_list = [file.name for file in file_list]
        byte_size = original_file[0].size
        bit_size = byte_size*8
        print(bit_size)
        dna_sequences_list = [self.get_dna_sequence(file) for file in file_list]
        user = request.user
        key = {
            "bit_size": bit_size,
            "file_name_list": file_name_list,
            "dna_sequences_list": dna_sequences_list,

        }
        try:
            if user.id is not None:
                key['user_id'] = user.id
                def task():
                    # 这里执行你的计算任务
                    print(f'Starting task at {datetime.now()} in thread: {current_thread().name}')
                    test02(**key)
                    # time.sleep(30)
                    print(f'Ending task at {datetime.now()} in thread: {current_thread().name}')
                thread_pool.submit(task)
                # thread.join()
                return JsonResponse({
                    "code": 200,
                    "data":"",
                    "msg": "msg",
                })
            else:
                data = test03(**key)
                return JsonResponse({
                    "code": 200,
                    "data": json.dumps(data),
                    "msg": "msg",
                })
        except Exception as e:
            print(e)
            return JsonResponse({
                "code": 400,
            })



# 查询是否有计算的结果记录
class GetData(View):
    def get(self, request):
        try:
            user_id = request.user.id
            key = request.GET['key']
            # user_id = 1

            tbResultModel = TbResult.objects.filter(user_id=user_id, dataType=key).order_by("-id").first()
            data = tbResultModel.result if tbResultModel.result else None
            print("##########", data)
            return JsonResponse({
                "code": 200,
                "data": json.loads(data),
                "msg": "OK"
            })
        except Exception as e:
            print(e)
            return JsonResponse({
                "code": 999999,
                "data": None,
                "msg": "None"
            })


class GetTask(View):
    def post(self, request):
        key_id = request.POST.get('id')
        user = request.user
        if key_id is not None and key_id != '' and user.id is not None:
            TbResult.objects.filter(id=key_id).update(is_valid=True)
        result = None
        if user.id is not None:
            mo = TbResult.objects.filter(is_valid=False, user_id=user.id).order_by("-id").values("id", "error", "dataType").all()
            result = []
            for i in mo:
                result.append({
                    "id": i.get("id"),
                    "message": i.get("dataType"),
                    "error": i.get("error"),
                    "dataType": i.get("dataType") if i.get("dataType") else None,
                })
        return JsonResponse({
            "code": 200,
            "data": result,
            "msg": "OK"
        })


class GetTaskDataInfo(View):
    def post(self, request):
        key_id = request.POST.get('id')
        mo = TbResult.objects.filter(id=key_id).first()
        result = {
                "id": mo.id,
                "message": mo.dataType,
                "error": mo.error,
                "dataType": mo.dataType if mo.dataType else None,
                "data": mo.result if mo.result else None
            }
        return JsonResponse({
            "code": 200,
            "data": result,
            "msg": "OK"
        })
