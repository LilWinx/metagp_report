from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib import utils, colors
from reportlab.platypus import Table, TableStyle, SimpleDocTemplate, Spacer, Image
from io import BytesIO

def pdf_template(outdir, in_dict):
    outpath = outdir + "/report.pdf"
    c = canvas.Canvas(outpath, pagesize=A4)

    # Letter head stuff
    logo_img = utils.ImageReader(in_dict.get("py_logo_ph", ""))
    metagp_img = utils.ImageReader(in_dict.get("py_metagplogo_ph", ""))

    c.drawImage(logo_img, 10, A4[1] - 54, width=110, height=44, mask="auto")
    c.drawImage(metagp_img, (A4[0] / 2) - 55, A4[1] - 54, width=110, height=30, mask="auto")

    address = [
        "Centre for Infectious Diseases and Microbiology Laboratory Services",
        "Level 3, NSW Health Pathology, ICPMR, Westmead Hospital",
        "Westmead, NSW, 2145, Australia",
    ]
    address_y = A4[1] - 30
    c.setFont("Helvetica", 5)
    for line in address:
        c.drawRightString(A4[0] - 10, address_y, line)
        address_y -= 6

    # Title
    title = "Metagenomics Report"
    c.setFont("Helvetica", 12)
    c.drawCentredString(A4[0] / 2, A4[1] - 80, title)
    c.translate(0, -60) # add space after
    
    # Patient Info table
    pt_metadata = [['Patient Name:', 'py_pn_ph', 'MRN:', 'py_mrn_ph'],
                    ["Accession:", "py_acc_ph", "Date of Collection:", "py_doc_ph"],
                    ["Date of Birth:", "py_dob_ph", "Request Date:", "py_reqdate_ph"],
                    ["Sequence Date:","py_seqdate_ph","Run Name:","py_runname_ph"],
                    ["Date of Report:","py_repdate_ph", "WGS ID:", "py_wgsid_ph"]]
    for row_idx, row in enumerate(pt_metadata):
        for col_idx, cell in enumerate(row):
            if cell in in_dict:
                pt_metadata[row_idx][col_idx] = in_dict[cell]

    c.setFont("Helvetica", 9)
    table = Table(pt_metadata, rowHeights=4*3.5)
    
    style = TableStyle([
        ('FONTNAME', (0, 0), (4, 4), 'Helvetica'),
        ('FONTSIZE', (0, 0), (4, 4), 8),
        ('RIGHTPADDING', (1, 0), (1, -1), 40),
        ('BOTTOMPADDING', (0, 0), (4, 4), 1),
    ])
    
    table.setStyle(style)
    table.wrapOn(c, 0,0)
    table.drawOn(c, A4[0] / 2 - 175, A4[1] - 110)
    
    # draw dividing line
    x_line, y_line = 50, A4[1] - 130
    x_line2, y_line2 = 550, A4[1] - 130
    c.line(x_line, y_line, x_line2, y_line2)

    # Results section  
    clinical_results_list = [
        "Date of Illness Onset: 'py_date_onset_ph'",
        "Suspected Pathogen: 'py_suspathogen_ph'",
        "Sample Type: 'py_sampletype_ph'",
        "Clinician Comments: 'py_clincomm_ph'",
    ]

    meta_summary_list = [
        "Interpretation: 'py_riinterpret_ph'",
        "Reason: 'py_rireason_ph'"
    ]
    line_count = 0
    header = "Clinical Request Information:"
    header2 = "Metagenomics Summary:"
    c.setFont("Helvetica-Bold", 8)
    c.drawString(50, A4[1]- 150, header)
    text_width = c.stringWidth(header, "Helvetica-Bold", 8)
    c.line(50, A4[1] - 152, 50 + text_width, A4[1] - 152)
    
    c.setFont("Helvetica", 8)

    # Your clinical results
    clinical_buffer = A4[1] - 170
    for result in clinical_results_list:
        placeholder = result.split(": ")[1].strip("'")
        if placeholder in in_dict:
            result = result.replace(f"'{placeholder}'", in_dict[placeholder])
            numlines = draw_wrapped_text(c, result, 50, clinical_buffer, 500)
            clinical_buffer -= 15
            line_count += numlines
    
    # Your meta summary
    summary_buffer = clinical_buffer - (line_count * 8)
    c.setFont("Helvetica-Bold", 8)
    c.drawString(50, summary_buffer, header2)
    text_width2 = c.stringWidth(header2, "Helvetica-Bold", 8)
    c.line(50, summary_buffer - 2, 50 + text_width2, summary_buffer - 2)

    summary_buffer = summary_buffer - 15
    for summary in meta_summary_list:
        placeholder = summary.split(": ")[1].strip("'")
        if placeholder in in_dict:
            summary = summary.replace(f"'{placeholder}'", in_dict[placeholder])
            draw_wrapped_text(c, summary, 50, summary_buffer, 500)
            summary_buffer -= 15
    
    c.showPage()
    c.save()

def draw_wrapped_text(c, text, x, y, width):
    # Set the font and font size
    font_name="Helvetica"
    font_size=8
    c.setFont(font_name, font_size)

    # Create a TextObject for more control
    text_object = c.beginText(x, y)
    text_object.setFont(font_name, font_size)

    # Split the text into lines manually to fit within the specified width
    lines = []
    current_line = ""
    words = text.split()
    
    for word in words:
        text_width = c.stringWidth(current_line + " " + word, font_name, font_size)
        if text_width < width:
            current_line += " " + word
        else:
            lines.append(current_line.strip())
            current_line = word
    
    # Add the last line
    lines.append(current_line.strip())

    # Draw each line
    for line in lines:
        text_object.textLine(line)

    # Draw the text on the canvas
    c.drawText(text_object)

    # Draw each line
    for line in lines:
        text_object.textLine(line)

    # Draw the text on the canvas
    c.drawText(text_object)
    return len(lines)