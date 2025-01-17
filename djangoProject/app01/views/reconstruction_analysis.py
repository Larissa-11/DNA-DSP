import json
import io
from threading import Thread
from django.views import View
from django.http import JsonResponse
from app01.Reconstruction_analysis.compute_performances import *
from app01.models import TbResult
from app01.views.algorithm_evaluation import thread_pool
from app01.views.algorithm_evaluation import notification_queue
from app01.views.errorInfoList import ERROR_INFO_LIST
from djangoProject.tools.webSocket import send_notification


class ClustringAnalysis(View):
    def post(self, request):
        try:
            true_file_list = request.FILES.getlist('true_file')
            pred_file_list = request.FILES.getlist('pred_file')

            user_id = request.user.id
            # user_id = 1  # 示例用户ID

            if not true_file_list or not pred_file_list:
                return JsonResponse({
                    "code": 400,
                    "msg": "No files provided",
                })

            # 将文件保存到内存
            true_files = self.handle_uploaded_files_in_memory(true_file_list)
            pred_files = self.handle_uploaded_files_in_memory(pred_file_list)

            def clustering_analysis_wrapper(true_files, pred_files, user_id):
                try:
                    # 调用实际的分析函数
                    ClustringAnalysisFunc(true_files, pred_files, user_id)
                except Exception as e:
                    print(f"Error in thread: {e}")

            data = {
                "true_files": true_files,
                "pred_files": pred_files,
                "user_id": user_id
            }
            # 启动子线程进行分析
            # thread = Thread(target=clustering_analysis_wrapper, args=(true_files, pred_files, user_id))
            # thread.start()
            if user_id:
                thread_pool.submit(clustering_analysis_wrapper, **data)
                return JsonResponse({
                    "code": 200,
                    "msg": "Success",
                })
            else:
                data = ClustringAnalysisFunc_01(true_files=true_files, pred_files=pred_files)
                return JsonResponse({
                    "code": 200,
                    "data": json.dumps(data),
                    "msg": "Success",
                })
        except Exception as e:
            return JsonResponse({
                "code": 500,
                "msg": str(e),
            })

    def handle_uploaded_files_in_memory(self, file_list):
        files_in_memory = []
        for file in file_list:
            file_in_memory = io.BytesIO(file.read())
            file_in_memory.name = file.name
            files_in_memory.append(file_in_memory)
        return files_in_memory


def ClustringAnalysisFunc(true_files, pred_files, user_id):
    try:
        # 获取第一个真实文件
        true_file = true_files[0]
        true_labels_new = getting_true_labels_new(true_file)
        results = []
        for pred_file in pred_files:
            pred_labels_new = getting_cluster_labels_new(pred_file)
            nmi_new = compute_NMI(true_labels_new, pred_labels_new)
            ami_new = compute_AMI(true_labels_new, pred_labels_new)
            ari_new = compute_ARI(true_labels_new, pred_labels_new)
            fmi_new = compute_FMI(true_labels_new, pred_labels_new)
            homo_new = compute_HOMO(true_labels_new, pred_labels_new)
            comp_new = compute_COMP(true_labels_new, pred_labels_new)

            results.append({
                "pred_file": pred_file.name,
                "NMI_new": round(nmi_new, 3),
                "AMI_new": round(ami_new, 3),
                "ARI_new": round(ari_new, 3),
                "FMI_new": round(fmi_new, 3),
                "HOMO_new": round(homo_new, 3),
                "COMP_new": round(comp_new, 3)
            })

        tmp_data = {
            "file_name": [],
            "NMI_new": [],
            "AMI_new": [],
            "ARI_new": [],
            "FMI_new": [],
            "HOMO_new": [],
            "COMP_new": [],
        }
        for i in results:
            tmp_data["file_name"].append(i["pred_file"])
            tmp_data["NMI_new"].append(i["NMI_new"])
            tmp_data["AMI_new"].append(i["AMI_new"])
            tmp_data["ARI_new"].append(i["ARI_new"])
            tmp_data["FMI_new"].append(i["FMI_new"])
            tmp_data["HOMO_new"].append(i["HOMO_new"])
            tmp_data["COMP_new"].append(i["COMP_new"])

        tmp_model = TbResult.objects.create(result=json.dumps(tmp_data), user_id=user_id, dataType="ClustringAnalysis")
        send_notification(tmp_model.id, "ClustringAnalysis", "ClustringAnalysis")
        return tmp_data
    except Exception as e:
        error_info = str(e)
        key_name = "ClustringAnalysisError"
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

def ClustringAnalysisFunc_01(true_files, pred_files):
    try:
        # 获取第一个真实文件
        true_file = true_files[0]
        true_labels_new = getting_true_labels_new(true_file)
        results = []
        for pred_file in pred_files:
            pred_labels_new = getting_cluster_labels_new(pred_file)
            nmi_new = compute_NMI(true_labels_new, pred_labels_new)
            ami_new = compute_AMI(true_labels_new, pred_labels_new)
            ari_new = compute_ARI(true_labels_new, pred_labels_new)
            fmi_new = compute_FMI(true_labels_new, pred_labels_new)
            homo_new = compute_HOMO(true_labels_new, pred_labels_new)
            comp_new = compute_COMP(true_labels_new, pred_labels_new)

            results.append({
                "pred_file": pred_file.name,
                "NMI_new": round(nmi_new, 3),
                "AMI_new": round(ami_new, 3),
                "ARI_new": round(ari_new, 3),
                "FMI_new": round(fmi_new, 3),
                "HOMO_new": round(homo_new, 3),
                "COMP_new": round(comp_new, 3)
            })

        tmp_data = {
            "file_name": [],
            "NMI_new": [],
            "AMI_new": [],
            "ARI_new": [],
            "FMI_new": [],
            "HOMO_new": [],
            "COMP_new": [],
        }
        for i in results:
            tmp_data["file_name"].append(i["pred_file"])
            tmp_data["NMI_new"].append(i["NMI_new"])
            tmp_data["AMI_new"].append(i["AMI_new"])
            tmp_data["ARI_new"].append(i["ARI_new"])
            tmp_data["FMI_new"].append(i["FMI_new"])
            tmp_data["HOMO_new"].append(i["HOMO_new"])
            tmp_data["COMP_new"].append(i["COMP_new"])
        return tmp_data
    except Exception as e:
        print("error:", e)

