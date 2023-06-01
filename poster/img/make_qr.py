import qrcode
from qrcode.image.styledpil import StyledPilImage

qrs_to_make = {
    'github-qr.png': 'https://github.com/jacopok/lgwa-skyloc',
    'arxiv-qr.png': 'https://arxiv.org/abs/2010.13726'
}

for filename, url in qrs_to_make.items():
    qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_L)
    qr.add_data(url)
    img = qr.make_image(image_factory=StyledPilImage)
    img.save(filename)