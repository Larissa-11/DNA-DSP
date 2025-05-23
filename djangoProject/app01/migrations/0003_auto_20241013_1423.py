# Generated by Django 3.2.25 on 2024-10-13 14:23

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('app01', '0002_auto_20240128_2129'),
    ]

    operations = [
        migrations.AddField(
            model_name='dnainfo',
            name='bits',
            field=models.IntegerField(default=0),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='dnainfo',
            name='content_type',
            field=models.CharField(max_length=255, verbose_name='文件类型'),
        ),
        migrations.CreateModel(
            name='TbResult',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('result', models.JSONField(blank=True, null=True, verbose_name='结果存储字段')),
                ('dataType', models.CharField(max_length=255, verbose_name='类型')),
                ('error', models.CharField(max_length=255, verbose_name='类型')),
                ('is_valid', models.BooleanField(default=False, verbose_name='逻辑删除字段，默认为0')),
                ('user', models.ForeignKey(null=True, on_delete=django.db.models.deletion.RESTRICT, to=settings.AUTH_USER_MODEL, verbose_name='用户外键。')),
            ],
            options={
                'db_table': 'tb_results',
            },
        ),
    ]
