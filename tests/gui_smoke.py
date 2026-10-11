"""需真实/虚拟显示环境。实际调用 GUI 生成与保存按钮，验证失败后清空旧图。"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import tkinter as tk
from tkinter import filedialog
from PIL import Image, ImageGrab
import chem_drawer

original = tk.Tk
errors = []
output = Path('output')
output.mkdir(exist_ok=True)

def factory():
    root = original()
    def descendants(widget):
        result = []
        for child in widget.winfo_children():
            result.append(child)
            result.extend(descendants(child))
        return result
    def start():
        try:
            widgets = descendants(root)
            if '--molecule' in sys.argv or '--straight' in sys.argv:
                entry = next(w for w in widgets if w.winfo_class() == 'TEntry')
                entry.delete(0, 'end')
                entry.insert(0, 'SMILES:C[N+](C)(C)CC(=O)NN.[Cl-]' if '--straight' in sys.argv else 'Nα-(2,4-dinitro-5-fluorophenyl)-L-alaninamide')
                next(w for w in widgets if w.winfo_class() == 'TRadiobutton' and w.cget('text') == '单个分子').invoke()
            generate = next(w for w in widgets if w.winfo_class() == 'TButton' and w.cget('text') == '生成图片')
            generate.invoke()
            root.after(100, lambda: finish(widgets, 0))
        except Exception as exc:
            errors.append(exc); root.destroy()
    def finish(widgets, attempts):
        try:
            save = next(w for w in widgets if w.winfo_class() == 'TButton' and w.cget('text') == '保存图片')
            if str(save.cget('state')) == 'disabled':
                assert attempts < 100, 'GUI 未完成生成'
                root.after(100, lambda: finish(widgets, attempts+1))
                return
            for fmt in ['PNG', 'SVG']:
                button = next(w for w in widgets if w.winfo_class() == 'TRadiobutton' and w.cget('text') == fmt)
                button.invoke()
                filedialog.asksaveasfilename = lambda **kw: str(output/f'gui-saved.{fmt.lower()}')
                save.invoke()
            Image.open(output/'gui-saved.png').load()
            assert (output/'gui-saved.svg').read_text().startswith('<svg')
            assert next(w for w in widgets if w.winfo_class() == 'Canvas').find_all()
            root.update()
            ImageGrab.grab(xdisplay=':99').save(output/'gui-screen.png')
            entry = next(w for w in widgets if w.winfo_class() == 'TEntry')
            entry.delete(0, 'end'); entry.insert(0, 'C2H6O' if '--molecule' in sys.argv or '--straight' in sys.argv else 'C2H6O -> C2H6O')
            next(w for w in widgets if w.winfo_class() == 'TButton' and w.cget('text') == '生成图片').invoke()
            root.after(300, lambda: rejected(widgets, save))
        except Exception as exc:
            errors.append(exc); root.destroy()
    def rejected(widgets, save):
        try:
            label = next(w for w in widgets if w.winfo_class() == 'TLabel' and w.cget('textvariable'))
            assert '不能唯一确定' in root.getvar(label.cget('textvariable'))
            assert str(save.cget('state')) == 'disabled'
            assert not next(w for w in widgets if w.winfo_class() == 'Canvas').find_all()
        except Exception as exc:
            errors.append(exc)
        root.destroy()
    root.after(300, start)
    return root

tk.Tk = factory
chem_drawer.launch_gui()
if errors:
    raise errors[0]
print('GUI: preview, PNG/SVG save, ambiguity error and stale-image clearing passed')
