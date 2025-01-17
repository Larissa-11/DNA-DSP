from django.contrib.auth.models import User
from django.db import models


class DnaInfo(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='dna_info', db_index=True)
    file_name = models.CharField(max_length=64, verbose_name="文件名")
    suffix_name = models.CharField(max_length=16, verbose_name="文件后缀")
    content_type = models.CharField(max_length=255, verbose_name="文件类型")
    bit_size = models.IntegerField(verbose_name="DNA大小")
    segment_length = models.IntegerField(verbose_name="片段长度", help_text="这是加入索引纠错的长度")
    index_length = models.IntegerField(verbose_name="索引长度")
    total_count = models.IntegerField(verbose_name="解码包数")
    bits = models.IntegerField()
    prime_segment_length = models.IntegerField(verbose_name="最初片段长度", help_text="编码页面输入的片段长度")

    class Meta:
        db_table = 'dna_info'
        verbose_name = 'DNA信息'
        constraints = [
            models.UniqueConstraint(fields=['user', 'file_name'], name='unique_user_file_name')
        ]


class TbResult(models.Model):
    result = models.JSONField(null=True, blank=True, verbose_name="结果存储字段")
    dataType = models.CharField(max_length=255, verbose_name="类型")
    error = models.CharField(max_length=255, verbose_name="类型")
    is_valid = models.BooleanField(default=False,verbose_name="逻辑删除字段，默认为0")
    user = models.ForeignKey(User, on_delete=models.RESTRICT, null=True, verbose_name="用户外键。")

    class Meta:
        db_table = 'tb_results'
