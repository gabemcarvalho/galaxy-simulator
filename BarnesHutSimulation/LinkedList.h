#pragma once

template<class TYPE>
class List
{
public:
	typedef List<TYPE> ContainerType;
	typedef TYPE* Type;

	template<typename ReturnType>
	class BaseIterator
	{
	public:
		BaseIterator(const ContainerType* container = 0, ReturnType curr = 0) : m_Container(container), m_Curr(curr) {}

		ReturnType	GetValue() const { return m_Curr; }
		bool		AtEnd() const { return m_Curr == m_Container->GetEnd(); }
		bool		AtStart() const { return m_Curr == m_Container->GetStart(); }
		void		Next() { m_Curr = m_Curr->next; }
		bool		Done() const { return m_Curr == 0 || m_Container == 0; }

	protected:
		const ContainerType* GetContainer() const { return m_Container; }

	private:
		const ContainerType* m_Container;
		ReturnType m_Curr;
	};

	typedef BaseIterator<Type> Iterator;

	List() : m_Head(0), m_Tail(0), m_iSize(0) { }
	~List() {}

	Iterator GetIterator() { return Iterator(this, m_Head); }

	void Clear() { m_Head = m_Tail = 0; m_iSize = 0; }

	TYPE* GetStart() { return m_Head; }
	void SetStart(TYPE* pEntry) { m_Head = pEntry; }
	TYPE* GetEnd() { return m_Tail; }
	void SetEnd(TYPE* pEntry) { m_Tail = pEntry; }
	int GetSize() { return m_iSize; }

	void AddEnd(TYPE* element)
	{
		element->next = 0;

		if (m_Tail != 0)
		{
			// put this node after the current tail
			m_Tail->next = element;
			m_Tail = element;
		}
		else if (m_Head != 0)
		{
			// put this node after the head
			m_Head->next = element;
			m_Tail = element;
		}
		else
		{
			m_Tail = element;
			m_Head = element;
		}

		m_iSize++;
	}

	TYPE* RemoveElement(TYPE* element)
	{
		// trying to remove from an empty list.
		if (m_Head == 0) return 0;
		
		m_iSize--;

		// special case for head removal (no previous)
		if (m_Head == element)
		{
			// if the tail is the head (only one element in the list), nullify the tail.
			if (m_Tail != 0 && m_Tail == m_Head)
			{
				m_Tail = 0;
			}

			m_Head = m_Head->next;

			return 0;
		}

		TYPE* current = m_Head->next;
		TYPE* prev = m_Head;

		while (current != 0)
		{
			if (current == element)
			{
				prev->next = current->next;

				if (m_Tail != 0 && m_Tail == current)
				{
					m_Tail = prev;
				}
				return prev;
			}
			prev = current;
			current = current->next;
		}
		return 0;
	}

private:
	TYPE* m_Head;
	TYPE* m_Tail;
	int m_iSize;
};