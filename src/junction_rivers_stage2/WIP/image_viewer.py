import os
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk

class ImageViewer(tk.Tk):
    def __init__(self, folder_path):
        super().__init__()
        self.folder_path = folder_path
        self.image_files = [f for f in os.listdir(folder_path) if f.endswith('.png')]
        self.image_index = 0

        self.title("Image Viewer")
        self.geometry("800x600")

        self.canvas = tk.Canvas(self, bg='white')
        self.canvas.pack(fill=tk.BOTH, expand=True)

        self.scroll_y = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.canvas.yview)
        self.scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        self.scroll_x = ttk.Scrollbar(self, orient=tk.HORIZONTAL, command=self.canvas.xview)
        self.scroll_x.pack(side=tk.BOTTOM, fill=tk.X)

        self.canvas.configure(yscrollcommand=self.scroll_y.set, xscrollcommand=self.scroll_x.set)
        self.canvas.bind('<Configure>', self.resize_image)
        self.bind('<space>', self.next_image)
        self.bind('<MouseWheel>', self.zoom_image)

        # Bind middle mouse button for panning
        self.canvas.bind('<ButtonPress-2>', self.start_pan)
        self.canvas.bind('<B2-Motion>', self.do_pan)

        self.load_image()

        self.pan_start_x = 0
        self.pan_start_y = 0

    def load_image(self):
        if not self.image_files:
            return

        image_path = os.path.join(self.folder_path, self.image_files[self.image_index])
        self.image = Image.open(image_path)
        self.tk_image = ImageTk.PhotoImage(self.image)
        self.canvas_image = self.canvas.create_image(0, 0, anchor='nw', image=self.tk_image)
        self.canvas.config(scrollregion=self.canvas.bbox(tk.ALL))

    def next_image(self, event):
        self.image_index = (self.image_index + 1) % len(self.image_files)
        self.load_image()

    def resize_image(self, event):
        self.canvas.config(scrollregion=self.canvas.bbox(tk.ALL))

    def zoom_image(self, event):
        if event.delta > 0:
            self.image = self.image.resize((int(self.image.width * 1.1), int(self.image.height * 1.1)), Image.ANTIALIAS)
        elif event.delta < 0:
            self.image = self.image.resize((int(self.image.width * 0.9), int(self.image.height * 0.9)), Image.ANTIALIAS)

        self.tk_image = ImageTk.PhotoImage(self.image)
        self.canvas.itemconfig(self.canvas_image, image=self.tk_image)
        self.canvas.config(scrollregion=self.canvas.bbox(tk.ALL))

    def start_pan(self, event):
        self.pan_start_x = event.x
        self.pan_start_y = event.y

    def do_pan(self, event):
        dx = event.x - self.pan_start_x
        dy = event.y - self.pan_start_y
        self.canvas.xview_scroll(-dx, 'units')
        self.canvas.yview_scroll(-dy, 'units')
        self.pan_start_x = event.x
        self.pan_start_y = event.y

if __name__ == "__main__":
    folder_path = "G:/Bungaban/PSCAD_Models/Bungaban_GW_GF1_PSCAD_Model_v0.14/results/2024_06_06_1458_results/CSR_S5255"
    app = ImageViewer(folder_path)
    app.mainloop()