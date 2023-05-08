import base64

def html_base64_encode(html):
    with open(html, "rb") as krona_html:
        krona_html_data = krona_html.read()
    
    krona_base64 = base64.b64encode(krona_html_data).decode("utf-8")
    return krona_base64