from channels.layers import get_channel_layer
from asgiref.sync import async_to_sync


def send_notification(key_id, message, type_tmp, data=None, error=None):
    channel_layer = get_channel_layer()
    async_to_sync(channel_layer.group_send)(
        'notifications',
        {
            "id": key_id,
            'type': 'send_notification',
            'message': message,
            'error': error,
            "dataType": type_tmp,
            "data": data
        }
    )
