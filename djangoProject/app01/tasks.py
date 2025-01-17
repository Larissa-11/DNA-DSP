# myapp/tasks.py
import json
import multiprocessing
import time
from concurrent.futures import ThreadPoolExecutor
from threading import Thread

from app01.encoding_analysis.GC_content_distribution import GcContentDistribution
from app01.encoding_analysis.MFE import Mfe
from app01.encoding_analysis.Net_hom_motifs import CalculateNetHomopolymerMotifs
from app01.models import TbResult
# from djangoProject.celery_config import app as celery_app
from celery import shared_task

from app01.views.algorithm_evaluation import notification_queue
from app01.views.errorInfoList import ERROR_INFO_LIST
from djangoProject.tools.webSocket import send_notification


# @shared_task
# def add(x, y):
#     for i in range(100000):
#         print(i)
#     return x + y


@shared_task
def test():
    try:
        for i in range(100):
            print(i, "#####")
        tmp_ = {
            "name": "Tom",
            "age": 30
        }
        TbResult.objects.create(result=json.dumps(tmp_), user_id=1)
    except Exception as e:
        print("error:", e)

# @shared_task
def test02(bit_size=None, dna_sequences_list=None, file_name_list=None, user_id=None):
    try:
        print("111111111111111111")
        # time.sleep(5)
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
        if user_id:
            tmp_model = TbResult.objects.create(result=json.dumps(data),
                                                user_id=user_id,
                                                dataType="Encode_Evalution")
            send_notification(tmp_model.id, "Encode_Evalution", "Encode_Evalution")
        return data
    except Exception as e:
        print("error:", e)
        error_info = str(e)
        key_name = "EncodeEvalutionError"
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

# @shared_task
def test03(bit_size=None, dna_sequences_list=None, file_name_list=None, user_id=None):
    try:
        print("111111111111111111")
        # time.sleep(5)
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

        return data
    except Exception as e:
        print("error:", e)

        return ""